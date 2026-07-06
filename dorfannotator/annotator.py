from __future__ import annotations

import csv
import gzip
import json
import logging
import re
import sqlite3
from pathlib import Path
from typing import Iterable

logger = logging.getLogger(__name__)

_SIMPLE_ALLELE = re.compile(r"^[ACGTacgt]+$")


def is_simple_allele(ref: str, alt: str) -> bool:
    return bool(_SIMPLE_ALLELE.match(ref)) and bool(_SIMPLE_ALLELE.match(alt))

from .models import DORF, Effect, OUTPUT_FIELDS, Transcript, Variant
from .sequence import (
    DEFAULT_START_CODONS,
    STOP_CODONS,
    apply_variant_to_sequence,
    kozak_context,
    kozak_strength,
    revcomp,
    chrom_aliases,
)


PRIORITY = {
    "dStart_lost": 0,
    "dStart_changed": 1,
    "dStop_lost": 2,
    "dStop_changed": 3,
    "dFrameshift": 4,
    "dStop_gained": 5,
    "dStart_gained": 6,
    "dInframe": 7,
    "dKozak_changed": 8,
    "dMissense": 9,
    "dSynonymous": 10,
}

KOZAK_RANK = {"strong": 0, "moderate": 1, "weak": 2}

CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def translate(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[index : index + 3], "X") for index in range(0, len(seq) - 2, 3))


def aa(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), "X")


class Annotator:
    def __init__(
        self,
        db_path: str | Path,
        *,
        start_codons: Iterable[str] = DEFAULT_START_CODONS,
        stop_codons: Iterable[str] = STOP_CODONS,
        bin_size: int = 10000,
        chromosomes: Iterable[str] | None = None,
        region: str | None = None,
        load_chromosomes: Iterable[str] | None = None,
        mane_only: bool = False,
        include_predicted: bool = False,
        evidence_only: bool = False,
    ) -> None:
        self.db_path = str(db_path)
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self.start_codons = {codon.upper() for codon in start_codons}
        self.stop_codons = {codon.upper() for codon in stop_codons}
        metadata = self._load_metadata()
        self.min_length = int(metadata.get("min_length", "30"))
        self.max_length = int(metadata.get("max_length", "303"))
        self.bin_size = int(metadata.get("bin_size", str(bin_size)))
        if metadata.get("start_codons"):
            self.start_codons = {codon.upper() for codon in metadata["start_codons"].split(",") if codon}
        if metadata.get("stop_codons"):
            self.stop_codons = {codon.upper() for codon in metadata["stop_codons"].split(",") if codon}
        self.region = parse_region(region) if region else None
        self.chromosomes = {alias for chrom in chromosomes or [] for alias in chrom_aliases(chrom)}
        if load_chromosomes is not None:
            load_source = list(load_chromosomes)
        else:
            load_source = list(chromosomes or [])
            if self.region:
                load_source.append(self.region[0])
        self._load_chromosomes = {alias for chrom in load_source for alias in chrom_aliases(chrom)}
        self.mane_only = mane_only
        self.include_predicted = include_predicted
        self.evidence_only = evidence_only
        self.transcripts = self._load_transcripts()
        self._transcript_bounds = self._build_transcript_bounds()
        self.dorfs = self._load_dorfs()
        self.dorfs_by_id = {dorf.dorf_id: dorf for dorf in self.dorfs}
        self.dorfs_by_transcript: dict[str, list[DORF]] = {}
        for dorf in self.dorfs:
            self.dorfs_by_transcript.setdefault(dorf.transcript_id, []).append(dorf)
        self._dorf_bin_index = self._load_bin_index("dorf_bins", "dorf_id")
        self._transcript_bin_index = self._load_bin_index("transcript_bins", "transcript_id")
        self._offset_cache: dict[str, dict[int, int]] = {}

    def _chrom_filter(self) -> tuple[str, tuple[str, ...]]:
        if not self._load_chromosomes:
            return "", ()
        placeholders = ",".join("?" for _ in self._load_chromosomes)
        return f" WHERE chrom IN ({placeholders})", tuple(self._load_chromosomes)

    def _load_bin_index(self, table: str, id_column: str) -> dict[tuple[str, int], list[str]]:
        clause, params = self._chrom_filter()
        index: dict[tuple[str, int], list[str]] = {}
        for row in self.conn.execute(f"SELECT chrom, bin, {id_column} FROM {table}{clause}", params):
            index.setdefault((row["chrom"], row["bin"]), []).append(row[id_column])
        return index

    def _pos_to_offset(self, transcript: Transcript) -> dict[int, int]:
        cached = self._offset_cache.get(transcript.transcript_id)
        if cached is None:
            positions = transcript.orf_genomic_positions or transcript.utr_genomic_positions
            cached = {pos: offset for offset, pos in enumerate(positions)}
            self._offset_cache[transcript.transcript_id] = cached
        return cached

    def _build_transcript_bounds(self) -> dict[str, tuple[int, int]]:
        bounds: dict[str, tuple[int, int]] = {}
        for transcript_id, transcript in self.transcripts.items():
            positions = transcript.orf_genomic_positions or transcript.utr_genomic_positions
            if positions:
                bounds[transcript_id] = (min(positions), max(positions))
        return bounds

    def _load_metadata(self) -> dict[str, str]:
        try:
            rows = self.conn.execute("SELECT key, value FROM metadata").fetchall()
        except sqlite3.OperationalError:
            return {}
        return {row["key"]: row["value"] for row in rows}

    def close(self) -> None:
        self.conn.close()

    def _load_transcripts(self) -> dict[str, Transcript]:
        clause, params = self._chrom_filter()
        rows = self.conn.execute(f"SELECT * FROM transcripts{clause}", params).fetchall()
        return {
            row["transcript_id"]: Transcript(
                transcript_id=row["transcript_id"],
                gene_id=row["gene_id"],
                gene_name=row["gene_name"],
                chrom=row["chrom"],
                strand=row["strand"],
                utr_sequence=row["utr_sequence"],
                utr_genomic_positions=tuple(json.loads(row["utr_genomic_positions"])),
                is_mane=bool(row["is_mane"]) if "is_mane" in row.keys() else False,
                orf_sequence=row["orf_sequence"] if "orf_sequence" in row.keys() else row["utr_sequence"],
                orf_genomic_positions=tuple(json.loads(row["orf_genomic_positions"])) if "orf_genomic_positions" in row.keys() else tuple(json.loads(row["utr_genomic_positions"])),
                utr_start_tx=row["utr_start_tx"] if "utr_start_tx" in row.keys() else 0,
                cds_tail_cds_start_offset=row["cds_tail_cds_start_offset"] if "cds_tail_cds_start_offset" in row.keys() else 0,
            )
            for row in rows
        }

    def _load_dorfs(self) -> list[DORF]:
        clause, params = self._chrom_filter()
        rows = self.conn.execute(f"SELECT * FROM dorfs{clause}", params).fetchall()
        dorfs = []
        for row in rows:
            dorfs.append(
                DORF(
                    dorf_id=row["dorf_id"],
                    transcript_id=row["transcript_id"],
                    gene_id=row["gene_id"],
                    gene_name=row["gene_name"],
                    chrom=row["chrom"],
                    strand=row["strand"],
                    start_tx=row["start_tx"],
                    end_tx=row["end_tx"],
                    genomic_start=row["genomic_start"],
                    genomic_end=row["genomic_end"],
                    start_codon=row["start_codon"],
                    stop_codon=row["stop_codon"],
                    length_nt=row["length_nt"],
                    kozak_sequence=row["kozak_sequence"],
                    kozak_strength=row["kozak_strength"],
                    orf_class=row["orf_class"] if "orf_class" in row.keys() else "dORF",
                    has_evidence=bool(row["has_evidence"]),
                    transcript_dorf_count_total=row["transcript_dorf_count_total"],
                    dorf_start_distance_from_cds_stop_nt=row["dorf_start_distance_from_cds_stop_nt"],
                    is_mane=bool(row["is_mane"]) if "is_mane" in row.keys() else False,
                )
            )
        return dorfs

    def annotate_batch(self, variants: Iterable[Variant]) -> list[Effect]:
        effects: list[Effect] = []
        for variant in variants:
            effects.extend(self.annotate_variant(variant))
        return effects

    def annotate_variant(self, variant: Variant) -> list[Effect]:
        if not self._variant_passes_filters(variant):
            return []
        candidates = self._candidate_dorfs(variant)
        effects: list[Effect] = []
        for dorf in candidates:
            if self.mane_only and not dorf.is_mane:
                continue
            transcript = self.transcripts[dorf.transcript_id]
            allele = self._variant_as_transcript_allele(transcript, variant)
            if allele is None:
                continue
            offset, ref, alt = allele
            effect = self._annotate_dorf(transcript, dorf, variant, offset, ref, alt)
            if effect is not None:
                effects.append(effect)
        gained = self._annotate_start_gained(variant)
        effects.extend(gained)
        effects.extend(self._annotate_stop_gained(variant))
        effects = self._collapse_predicted_stop_gained(effects)
        effects = [effect for effect in effects if self._effect_passes_output_mode(effect)]
        effects.sort(key=lambda item: (item.transcript_id, PRIORITY.get(item.consequence, 99), item.dorf_id))
        return effects

    def _effect_passes_output_mode(self, effect: Effect) -> bool:
        if self.include_predicted:
            return True
        if effect.has_evidence:
            return True
        if self.evidence_only:
            return False
        if effect.consequence == "dStart_gained":
            return effect.kozak_strength_alt == "strong"
        if effect.consequence == "dStop_gained":
            return (effect.kozak_strength_alt or effect.kozak_strength_ref) == "strong"
        if effect.consequence == "dKozak_changed":
            return effect.kozak_strength_alt == "strong"
        return False

    def _variant_passes_filters(self, variant: Variant) -> bool:
        if self.chromosomes and not any(alias in self.chromosomes for alias in chrom_aliases(variant.chrom)):
            return False
        if self.region:
            region_chrom, region_start, region_end = self.region
            if region_chrom not in chrom_aliases(variant.chrom):
                return False
            variant_end = variant.pos + max(1, len(variant.ref)) - 1
            if variant.pos > region_end or variant_end < region_start:
                return False
        return True

    def _bin_ids(self, variant: Variant) -> range:
        span_end = variant.pos + max(1, len(variant.ref)) - 1
        return range(variant.pos // self.bin_size - 1, span_end // self.bin_size + 2)

    def _candidate_dorfs(self, variant: Variant) -> list[DORF]:
        span_end = variant.pos + max(1, len(variant.ref)) - 1
        seen: set[str] = set()
        candidates: list[DORF] = []
        for alias in chrom_aliases(variant.chrom):
            for bin_id in self._bin_ids(variant):
                for dorf_id in self._dorf_bin_index.get((alias, bin_id), ()):
                    if dorf_id in seen:
                        continue
                    seen.add(dorf_id)
                    dorf = self.dorfs_by_id[dorf_id]
                    low = dorf.genomic_start - (6 if dorf.strand == "+" else 0)
                    high = dorf.genomic_end + (6 if dorf.strand == "-" else 0)
                    if low <= span_end and variant.pos <= high:
                        candidates.append(dorf)
        return candidates

    def _variant_as_transcript_allele(self, transcript: Transcript, variant: Variant) -> tuple[int, str, str] | None:
        pos_to_offset = self._pos_to_offset(transcript)
        ref = variant.ref.upper()
        alt = variant.alt.upper()
        seq = transcript.orf_sequence or transcript.utr_sequence
        minus = transcript.strand == "-"
        exonic_offsets: list[int] = []
        exonic_ref_indices: list[int] = []
        for index, base in enumerate(ref):
            offset = pos_to_offset.get(variant.pos + index)
            if offset is None:
                continue
            expected = revcomp(base) if minus else base
            if offset >= len(seq) or seq[offset].upper() != expected:
                return None
            exonic_offsets.append(offset)
            exonic_ref_indices.append(index)
        if not exonic_offsets:
            return None
        start_off, end_off = min(exonic_offsets), max(exonic_offsets)
        effective_ref = seq[start_off : end_off + 1]
        anchor_is_exonic = 0 in exonic_ref_indices
        if anchor_is_exonic:
            effective_alt = revcomp(alt) if minus else alt
        elif len(alt) == 1 and len(ref) > len(alt) and alt == ref[0]:
            effective_alt = ""
        else:
            return None
        return start_off, effective_ref, effective_alt


    def _genomic_bounds(self, transcript: Transcript, start_tx: int, end_tx: int) -> tuple[int, int]:
        all_positions = transcript.orf_genomic_positions or transcript.utr_genomic_positions
        positions = all_positions[max(0, start_tx) : min(end_tx, len(all_positions))]
        if not positions:
            return 0, 0
        return min(positions), max(positions)

    def _codon_start_pos(self, transcript: Transcript, codon_tx: int) -> int | None:
        positions = transcript.orf_genomic_positions or transcript.utr_genomic_positions
        if 0 <= codon_tx < len(positions):
            return positions[codon_tx]
        return None

    def _classify_orf(self, transcript: Transcript, start_tx: int, end_tx: int, utr_start_tx: int | None = None) -> str:
        utr_start = transcript.utr_start_tx if utr_start_tx is None else utr_start_tx
        stop_tx = end_tx - 3
        if start_tx >= utr_start and stop_tx >= utr_start:
            return "dORF"
        if start_tx < utr_start <= stop_tx:
            if (transcript.cds_tail_cds_start_offset + start_tx) % 3 != 0:
                return "doORF"
        return ""

    def _alt_utr_start(self, transcript: Transcript, offset: int, delta: int) -> int:
        return transcript.utr_start_tx + (delta if offset < transcript.utr_start_tx else 0)

    def _set_orf_output(self, effect: Effect, transcript: Transcript, start_tx: int, end_tx: int, utr_start_tx: int | None = None) -> None:
        effect.dorf_start, effect.dorf_end = self._genomic_bounds(transcript, start_tx, end_tx)
        utr_start = transcript.utr_start_tx if utr_start_tx is None else utr_start_tx
        effect.dist_cds = start_tx - utr_start
        alt_class = self._classify_orf(transcript, start_tx, end_tx, utr_start)
        if alt_class:
            effect.orf_class = alt_class

    def _alt_genomic_bounds(
        self,
        transcript: Transcript,
        offset: int,
        delta: int,
        alt_len: int,
        alt_start_tx: int,
        alt_end_tx: int,
    ) -> tuple[int, int]:
        positions = []
        for alt_index in (alt_start_tx, alt_end_tx - 1):
            pos = self._alt_index_to_genomic(transcript, offset, delta, alt_len, alt_index)
            if pos is not None:
                positions.append(pos)
        if positions:
            return min(positions), max(positions)
        return self._genomic_bounds(transcript, alt_start_tx, alt_end_tx)

    def _set_alt_orf_output(
        self,
        effect: Effect,
        transcript: Transcript,
        offset: int,
        delta: int,
        alt_len: int,
        start_tx: int,
        end_tx: int,
        utr_start_tx: int | None = None,
    ) -> None:
        effect.dorf_start, effect.dorf_end = self._alt_genomic_bounds(transcript, offset, delta, alt_len, start_tx, end_tx)
        utr_start = transcript.utr_start_tx if utr_start_tx is None else utr_start_tx
        effect.dist_cds = start_tx - utr_start
        alt_class = self._classify_orf(transcript, start_tx, end_tx, utr_start)
        if alt_class:
            effect.orf_class = alt_class

    def _class_change_suffix(self, effect: Effect, ref_class: str) -> str:
        return f";class={ref_class}>{effect.orf_class}" if effect.orf_class and effect.orf_class != ref_class else ""

    def _first_changed_codon(self, ref_orf: str, alt_orf: str) -> tuple[str, str, str, str]:
        limit = min(len(ref_orf), len(alt_orf)) - 2
        for index in range(0, limit, 3):
            ref_codon = ref_orf[index : index + 3]
            alt_codon = alt_orf[index : index + 3]
            if ref_codon != alt_codon:
                return ref_codon, alt_codon, aa(ref_codon), aa(alt_codon)
        return "", "", "", ""

    def _kozak_detail(
        self,
        seq: str,
        alt_seq: str,
        ref_start_tx: int,
        alt_start_tx: int,
        ref_strength: str,
        alt_strength: str,
    ) -> str:
        positions: list[str] = []
        changes: list[str] = []
        checks = (("-3", -3), ("+4", 3))
        for label, relative_index in checks:
            ref_index = ref_start_tx + relative_index
            alt_index = alt_start_tx + relative_index
            ref_base = seq[ref_index] if 0 <= ref_index < len(seq) else "N"
            alt_base = alt_seq[alt_index] if 0 <= alt_index < len(alt_seq) else "N"
            if ref_base != alt_base:
                positions.append(label)
                changes.append(f"{ref_base}>{alt_base}")
        pos_text = ",".join(positions) if positions else ""
        change_text = ",".join(changes) if changes else ""
        return f"kozak={ref_strength}>{alt_strength};pos={pos_text};change={change_text}"

    def _detail_with_kozak(self, detail: str, strength: str) -> str:
        if not strength or any(part.startswith("kozak=") for part in detail.split(";")):
            return detail
        if ";class=" in detail:
            before_class, class_suffix = detail.split(";class=", 1)
            return f"{before_class};kozak={strength};class={class_suffix}"
        return f"{detail};kozak={strength}"

    def _has_kozak_strength_change(
        self,
        seq: str,
        alt_seq: str,
        ref_start_tx: int,
        alt_start_tx: int,
        ref_strength: str,
        alt_strength: str,
    ) -> bool:
        if not ref_strength or not alt_strength or ref_strength == alt_strength:
            return False
        for relative_index in (-3, 3):
            ref_index = ref_start_tx + relative_index
            alt_index = alt_start_tx + relative_index
            ref_base = seq[ref_index] if 0 <= ref_index < len(seq) else "N"
            alt_base = alt_seq[alt_index] if 0 <= alt_index < len(alt_seq) else "N"
            if ref_base != alt_base:
                return True
        return False

    def _candidate_orf_rank(self, effect: Effect) -> tuple[int, int, int, int, int]:
        strength = effect.kozak_strength_alt or effect.kozak_strength_ref
        start_codon = effect.start_codon_alt or effect.start_codon_ref
        return (
            KOZAK_RANK.get(strength, 99),
            0 if start_codon == "ATG" else 1,
            effect.dorf_start_distance_from_cds_stop_nt,
            effect.dorf_start,
            effect.dorf_end,
        )

    def _best_candidate(self, effects: list[Effect]) -> Effect | None:
        return min(effects, key=self._candidate_orf_rank) if effects else None

    def _collapse_predicted_stop_gained(self, effects: list[Effect]) -> list[Effect]:
        grouped: dict[tuple[str, int | None, str], list[Effect]] = {}
        kept: list[Effect] = []
        for effect in effects:
            if effect.consequence == "dStop_gained" and not effect.has_evidence:
                grouped.setdefault(
                    (effect.transcript_id, effect.alternative_stop_position, effect.alternative_stop_codon),
                    [],
                ).append(effect)
            else:
                kept.append(effect)
        for group in grouped.values():
            best = self._best_candidate(group)
            if best is not None:
                kept.append(best)
        return kept

    def _base_effect(self, dorf: DORF, variant: Variant) -> Effect:
        return Effect(
            chrom=variant.chrom,
            pos=variant.pos,
            ref=variant.ref,
            alt=variant.alt,
            gene_id=dorf.gene_id,
            gene_name=dorf.gene_name,
            transcript_id=dorf.transcript_id,
            strand=dorf.strand,
            dorf_id=dorf.dorf_id,
            orf_class=dorf.orf_class,
            is_mane=dorf.is_mane,
            consequence="",
            transcript_dorf_count_total=dorf.transcript_dorf_count_total,
            dorf_start_distance_from_cds_stop_nt=dorf.dorf_start_distance_from_cds_stop_nt,
            start_codon_ref=dorf.start_codon,
            stop_codon_ref=dorf.stop_codon,
            kozak_strength_ref=dorf.kozak_strength,
            ref_dorf_length_nt=dorf.length_nt,
            has_evidence=dorf.has_evidence,
            dorf_start=dorf.genomic_start,
            dorf_end=dorf.genomic_end,
            dist_cds=dorf.dorf_start_distance_from_cds_stop_nt,
        )

    def _annotate_dorf(
        self, transcript: Transcript, dorf: DORF, variant: Variant, offset: int, ref: str, alt: str
    ) -> Effect | None:
        seq = transcript.orf_sequence or transcript.utr_sequence
        if seq[offset : offset + len(ref)].upper() != ref.upper():
            return None
        alt_seq = apply_variant_to_sequence(seq, offset, ref, alt)
        effect = self._base_effect(dorf, variant)
        delta = len(alt) - len(ref)
        alt_utr_start = self._alt_utr_start(transcript, offset, delta)
        var_start = offset
        var_end = offset + len(ref)
        stop_start = dorf.end_tx - 3
        overlaps_dorf = var_start < dorf.end_tx and var_end > dorf.start_tx
        overlaps_start = var_start < dorf.start_tx + 3 and var_end > dorf.start_tx
        overlaps_stop = var_start < dorf.end_tx and var_end > stop_start

        if len(alt) < len(ref) and var_start <= dorf.start_tx and var_end >= dorf.end_tx:
            effect.consequence = "dStart_lost"
            effect.start_codon_alt = ""
            effect.alt_dorf_length_nt = 0
            effect.detail = self._detail_with_kozak(f"start={dorf.start_codon}>;alt_start=false", effect.kozak_strength_ref)
            return effect

        start_shift = delta if var_end <= dorf.start_tx else 0
        alt_start = dorf.start_tx + start_shift
        alt_stop_start = stop_start + (delta if var_end <= stop_start else 0)
        alt_end = dorf.end_tx + (delta if var_end <= dorf.end_tx else 0)
        alt_start_codon = alt_seq[alt_start : alt_start + 3] if 0 <= alt_start <= len(alt_seq) - 3 else ""
        alt_stop_codon = alt_seq[alt_stop_start : alt_stop_start + 3] if 0 <= alt_stop_start <= len(alt_seq) - 3 else ""
        effect.start_codon_alt = alt_start_codon
        effect.stop_codon_alt = alt_stop_codon
        effect.kozak_strength_alt = kozak_strength(alt_seq, alt_start) if alt_start_codon else ""
        effect.alt_dorf_length_nt = alt_end - alt_start
        effect.dorf_length_delta_nt = effect.alt_dorf_length_nt - dorf.length_nt
        kozak_changed = self._has_kozak_strength_change(
            seq, alt_seq, dorf.start_tx, alt_start, effect.kozak_strength_ref, effect.kozak_strength_alt
        )

        if overlaps_start:
            if alt_start_codon not in self.start_codons:
                effect.consequence = "dStart_lost"
                self._fill_alternative_start(effect, transcript, dorf, alt_seq, alt_start, offset, delta, len(alt), alt_stop_start)
                detail = f"start={dorf.start_codon}>{alt_start_codon};alt_start={'true' if effect.alternative_start_found else 'false'}"
                if effect.alternative_start_found:
                    class_suffix = self._class_change_suffix(effect, dorf.orf_class)
                    detail += (
                        f";alt_start_codon={effect.alternative_start_codon}"
                        f";alt_start_pos={effect.alternative_start_position}"
                        f";alt_start_kozak={effect.alternative_start_kozak_strength}"
                        f";new_len={effect.alt_dorf_length_nt}"
                        f"{class_suffix}"
                    )
                effect.detail = self._detail_with_kozak(detail, effect.kozak_strength_ref)
                return effect
            if alt_start_codon != dorf.start_codon:
                effect.consequence = "dStart_changed"
                effect.detail = self._detail_with_kozak(
                    f"start={dorf.start_codon}>{alt_start_codon}",
                    effect.kozak_strength_alt or effect.kozak_strength_ref,
                )
                return effect

        if overlaps_stop:
            if alt_stop_codon not in self.stop_codons:
                effect.consequence = "dStop_lost"
                self._fill_alternative_stop(effect, transcript, dorf, alt_seq, alt_stop_start, offset, delta, len(alt))
                detail = f"stop={dorf.stop_codon}>{alt_stop_codon};alt_stop={'true' if effect.alternative_stop_found else 'false'}"
                if effect.alternative_stop_found:
                    class_suffix = self._class_change_suffix(effect, dorf.orf_class)
                    detail += (
                        f";alt_stop_codon={effect.alternative_stop_codon}"
                        f";alt_stop_pos={effect.alternative_stop_position}"
                        f";new_len={effect.alt_dorf_length_nt}"
                        f";extension={effect.extension_length_nt}"
                        f"{class_suffix}"
                    )
                effect.detail = self._detail_with_kozak(detail, effect.kozak_strength_ref)
                return effect
            if alt_stop_codon != dorf.stop_codon:
                effect.consequence = "dStop_changed"
                effect.detail = self._detail_with_kozak(f"stop={dorf.stop_codon}>{alt_stop_codon}", effect.kozak_strength_ref)
                return effect

        if overlaps_dorf and len(ref) != len(alt):
            if delta % 3 == 0:
                ref_orf = seq[dorf.start_tx : dorf.end_tx]
                alt_orf = alt_seq[alt_start:alt_end]
                new_stop = self._first_internal_new_stop(ref_orf, alt_orf)
                if new_stop is not None:
                    stop_tx = alt_start + new_stop
                    new_len = new_stop + 3
                    effect.consequence = "dStop_gained"
                    effect.alternative_stop_codon = alt_seq[stop_tx : stop_tx + 3]
                    effect.alternative_stop_position = self._alt_index_to_genomic(transcript, offset, delta, len(alt), stop_tx)
                    effect.alt_dorf_length_nt = new_len
                    effect.truncation_length_nt = dorf.length_nt - new_len
                    self._set_alt_orf_output(effect, transcript, offset, delta, len(alt), alt_start, stop_tx + 3, alt_utr_start)
                    effect.detail = (
                        f"new_stop={effect.alternative_stop_codon};new_stop_pos={effect.alternative_stop_position}"
                        f";new_len={new_len};truncation={effect.truncation_length_nt}{self._class_change_suffix(effect, dorf.orf_class)}"
                    )
                    effect.detail = self._detail_with_kozak(effect.detail, effect.kozak_strength_ref)
                    return effect
                effect.consequence = "dInframe"
                self._set_alt_orf_output(effect, transcript, offset, delta, len(alt), alt_start, alt_end, alt_utr_start)
                effect.detail = self._detail_with_kozak(
                    f"len={dorf.length_nt}>{effect.alt_dorf_length_nt}{self._class_change_suffix(effect, dorf.orf_class)}",
                    effect.kozak_strength_ref,
                )
                return effect
            effect.consequence = "dFrameshift"
            self._fill_first_downstream_stop(effect, transcript, dorf, alt_seq, dorf.start_tx, offset, delta, len(alt))
            if effect.alternative_stop_found:
                class_suffix = self._class_change_suffix(effect, dorf.orf_class)
                if effect.alt_dorf_length_nt < dorf.length_nt:
                    effect.truncation_length_nt = dorf.length_nt - effect.alt_dorf_length_nt
                    effect.detail = (
                        f"truncated=true;new_stop={effect.alternative_stop_codon}"
                        f";new_stop_pos={effect.alternative_stop_position};new_len={effect.alt_dorf_length_nt}"
                        f"{class_suffix}"
                    )
                    effect.detail = self._detail_with_kozak(effect.detail, effect.kozak_strength_ref)
                else:
                    effect.extension_length_nt = effect.alt_dorf_length_nt - dorf.length_nt
                    effect.detail = (
                        f"truncated=false;downstream_stop=true;new_stop={effect.alternative_stop_codon}"
                        f";new_stop_pos={effect.alternative_stop_position};new_len={effect.alt_dorf_length_nt}"
                        f";extension={effect.extension_length_nt}{class_suffix}"
                    )
                    effect.detail = self._detail_with_kozak(effect.detail, effect.kozak_strength_ref)
            else:
                effect.detail = self._detail_with_kozak("truncated=false;downstream_stop=false", effect.kozak_strength_ref)
            return effect

        if overlaps_dorf:
            ref_orf = seq[dorf.start_tx : dorf.end_tx]
            alt_orf = alt_seq[dorf.start_tx : dorf.end_tx]
            new_stop = self._first_internal_new_stop(ref_orf, alt_orf)
            if new_stop is not None:
                stop_tx = dorf.start_tx + new_stop
                new_len = new_stop + 3
                effect.consequence = "dStop_gained"
                effect.alternative_stop_codon = alt_seq[stop_tx : stop_tx + 3]
                effect.alternative_stop_position = self._codon_start_pos(transcript, stop_tx)
                effect.alt_dorf_length_nt = new_len
                effect.truncation_length_nt = dorf.length_nt - new_len
                self._set_orf_output(effect, transcript, dorf.start_tx, stop_tx + 3, alt_utr_start)
                effect.detail = (
                    f"new_stop={effect.alternative_stop_codon};new_stop_pos={effect.alternative_stop_position}"
                    f";new_len={new_len};truncation={effect.truncation_length_nt}{self._class_change_suffix(effect, dorf.orf_class)}"
                )
                effect.detail = self._detail_with_kozak(effect.detail, effect.kozak_strength_ref)
                return effect
            ref_codon, alt_codon, ref_aa, alt_aa = self._first_changed_codon(ref_orf, alt_orf)
            if kozak_changed:
                effect.consequence = "dKozak_changed"
                effect.detail = self._kozak_detail(
                    seq, alt_seq, dorf.start_tx, alt_start, effect.kozak_strength_ref, effect.kozak_strength_alt
                )
                return effect
            if translate(ref_orf) == translate(alt_orf):
                effect.consequence = "dSynonymous"
                effect.detail = self._detail_with_kozak(f"codon={ref_codon}>{alt_codon}", effect.kozak_strength_ref)
            else:
                effect.consequence = "dMissense"
                effect.detail = self._detail_with_kozak(f"codon={ref_codon}>{alt_codon};aa={ref_aa}>{alt_aa}", effect.kozak_strength_ref)
            return effect

        if kozak_changed:
            effect.consequence = "dKozak_changed"
            effect.detail = self._kozak_detail(
                seq, alt_seq, dorf.start_tx, alt_start, effect.kozak_strength_ref, effect.kozak_strength_alt
            )
            return effect
        return None

    def _alt_index_to_genomic(
        self, transcript: Transcript, offset: int, delta: int, alt_len: int, alt_index: int
    ) -> int | None:
        positions = transcript.orf_genomic_positions or transcript.utr_genomic_positions
        if alt_index < offset:
            ref_idx = alt_index
        elif alt_index < offset + alt_len:
            ref_idx = min(offset, len(positions) - 1)
        else:
            ref_idx = alt_index - delta
        if 0 <= ref_idx < len(positions):
            return positions[ref_idx]
        return None

    def _fill_alternative_start(
        self,
        effect: Effect,
        transcript: Transcript,
        dorf: DORF,
        seq: str,
        original_start: int,
        offset: int,
        delta: int,
        alt_len: int,
        alt_stop_start: int | None = None,
    ) -> None:
        stop_start = alt_stop_start if alt_stop_start is not None else dorf.end_tx - 3
        alt_end_tx = stop_start + 3
        matches: list[int] = []
        for start in range(original_start + 3, stop_start, 3):
            codon = seq[start : start + 3]
            if codon in self.start_codons and self._no_stop_between(seq, start + 3, stop_start):
                matches.append(start)
        if not matches:
            return
        start = min(matches, key=lambda item: abs(item - original_start))
        effect.alternative_start_found = True
        effect.alternative_start_codon = seq[start : start + 3]
        effect.alternative_start_position = self._alt_index_to_genomic(transcript, offset, delta, alt_len, start)
        effect.alternative_start_kozak_strength = kozak_strength(seq, start)
        effect.alt_dorf_length_nt = alt_end_tx - start
        effect.dorf_length_delta_nt = effect.alt_dorf_length_nt - dorf.length_nt
        self._set_alt_orf_output(effect, transcript, offset, delta, alt_len, start, alt_end_tx, self._alt_utr_start(transcript, offset, delta))

    def _fill_alternative_stop(
        self, effect: Effect, transcript: Transcript, dorf: DORF, seq: str, original_stop: int, offset: int, delta: int, alt_len: int
    ) -> None:
        for stop in range(original_stop + 3, len(seq) - 2, 3):
            if seq[stop : stop + 3] in self.stop_codons:
                effect.alternative_stop_found = True
                effect.alternative_stop_codon = seq[stop : stop + 3]
                effect.alternative_stop_position = self._alt_index_to_genomic(transcript, offset, delta, alt_len, stop)
                effect.extension_length_nt = stop - original_stop
                effect.alt_dorf_length_nt = stop + 3 - dorf.start_tx
                effect.dorf_length_delta_nt = effect.alt_dorf_length_nt - dorf.length_nt
                self._set_alt_orf_output(effect, transcript, offset, delta, alt_len, dorf.start_tx, stop + 3, self._alt_utr_start(transcript, offset, delta))
                return

    def _fill_first_downstream_stop(
        self, effect: Effect, transcript: Transcript, dorf: DORF, seq: str, start: int, offset: int, delta: int, alt_len: int
    ) -> None:
        for stop in range(start + 3, len(seq) - 2, 3):
            if seq[stop : stop + 3] in self.stop_codons:
                effect.alternative_stop_found = True
                effect.alternative_stop_codon = seq[stop : stop + 3]
                effect.alternative_stop_position = self._alt_index_to_genomic(transcript, offset, delta, alt_len, stop)
                effect.alt_dorf_length_nt = stop + 3 - start
                effect.dorf_length_delta_nt = effect.alt_dorf_length_nt - dorf.length_nt
                self._set_alt_orf_output(effect, transcript, offset, delta, alt_len, start, stop + 3, self._alt_utr_start(transcript, offset, delta))
                return

    def _no_stop_between(self, seq: str, start: int, stop: int) -> bool:
        return all(seq[index : index + 3] not in self.stop_codons for index in range(start, stop, 3))

    def _first_internal_new_stop(self, ref_orf: str, alt_orf: str) -> int | None:
        for index in range(3, len(alt_orf) - 3, 3):
            ref_codon = ref_orf[index : index + 3] if index + 3 <= len(ref_orf) else ""
            if ref_codon not in self.stop_codons and alt_orf[index : index + 3] in self.stop_codons:
                return index
        return None

    def _annotate_start_gained(self, variant: Variant) -> list[Effect]:
        effects: list[Effect] = []
        for transcript in self._candidate_transcripts(variant):
            transcript_effects: list[Effect] = []
            allele = self._variant_as_transcript_allele(transcript, variant)
            if allele is None:
                continue
            offset, ref, alt = allele
            seq = transcript.orf_sequence or transcript.utr_sequence
            if seq[offset : offset + len(ref)].upper() != ref.upper():
                continue
            alt_seq = apply_variant_to_sequence(seq, offset, ref, alt)
            delta = len(alt) - len(ref)
            alt_utr_start = self._alt_utr_start(transcript, offset, delta)
            for start in range(max(0, offset - 2), min(len(alt_seq) - 2, offset + len(alt)) + 1):
                if start < offset:
                    ref_idx = start
                elif start < offset + len(alt):
                    ref_idx = -1
                else:
                    ref_idx = start - delta
                ref_codon = seq[ref_idx : ref_idx + 3] if 0 <= ref_idx and ref_idx + 3 <= len(seq) else ""
                alt_codon = alt_seq[start : start + 3]
                if ref_codon in self.start_codons or alt_codon not in self.start_codons:
                    continue
                if any(d.start_tx == start for d in self.dorfs_by_transcript.get(transcript.transcript_id, [])):
                    continue
                stop = self._next_stop(alt_seq, start + 3)
                if stop is None:
                    continue
                length = stop + 3 - start
                if length < self.min_length or length > self.max_length:
                    continue
                orf_class = self._classify_orf(transcript, start, stop + 3, alt_utr_start)
                if not orf_class:
                    continue
                bounds = self._alt_genomic_bounds(transcript, offset, delta, len(alt), start, stop + 3)
                if bounds == (0, 0):
                    continue
                dummy = DORF(
                    dorf_id=f"{transcript.transcript_id}:{orf_class}:novel_start:{variant.pos}:{variant.alt}",
                    transcript_id=transcript.transcript_id,
                    gene_id=transcript.gene_id,
                    gene_name=transcript.gene_name,
                    chrom=transcript.chrom,
                    strand=transcript.strand,
                    start_tx=start,
                    end_tx=stop + 3,
                    genomic_start=bounds[0],
                    genomic_end=bounds[1],
                    start_codon=alt_codon,
                    stop_codon=alt_seq[stop : stop + 3],
                    length_nt=length,
                    kozak_sequence=kozak_context(alt_seq, start),
                    kozak_strength=kozak_strength(alt_seq, start),
                    orf_class=orf_class,
                    has_evidence=False,
                    transcript_dorf_count_total=self._evidence_orf_count(transcript.transcript_id),
                    dorf_start_distance_from_cds_stop_nt=start - alt_utr_start,
                    is_mane=transcript.is_mane,
                )
                effect = self._base_effect(dummy, variant)
                effect.consequence = "dStart_gained"
                effect.start_codon_ref = ref_codon
                effect.start_codon_alt = alt_codon
                effect.stop_codon_alt = dummy.stop_codon
                effect.kozak_strength_alt = dummy.kozak_strength
                effect.has_evidence = False
                effect.alt_dorf_length_nt = dummy.length_nt
                stop_pos = self._alt_index_to_genomic(transcript, offset, delta, len(alt), stop)
                effect.detail = f"new_start={alt_codon};new_stop={dummy.stop_codon};new_stop_pos={stop_pos};new_len={length};kozak={dummy.kozak_strength}"
                transcript_effects.append(effect)
            best = self._best_candidate(transcript_effects)
            if best is not None:
                effects.append(best)
        return effects

    def _annotate_stop_gained(self, variant: Variant) -> list[Effect]:
        effects: list[Effect] = []
        for transcript in self._candidate_transcripts(variant):
            allele = self._variant_as_transcript_allele(transcript, variant)
            if allele is None:
                continue
            offset, ref, alt = allele
            seq = transcript.orf_sequence or transcript.utr_sequence
            if seq[offset : offset + len(ref)].upper() != ref.upper():
                continue
            alt_seq = apply_variant_to_sequence(seq, offset, ref, alt)
            delta = len(alt) - len(ref)
            alt_utr_start = self._alt_utr_start(transcript, offset, delta)
            for stop in range(max(0, offset - 2), min(len(alt_seq) - 2, offset + len(alt)) + 1):
                if stop < offset:
                    ref_idx = stop
                elif stop < offset + len(alt):
                    ref_idx = -1
                else:
                    ref_idx = stop - delta
                ref_codon = seq[ref_idx : ref_idx + 3] if 0 <= ref_idx and ref_idx + 3 <= len(seq) else ""
                alt_codon = alt_seq[stop : stop + 3]
                if ref_codon in self.stop_codons or alt_codon not in self.stop_codons:
                    continue
                if any(
                    d.start_tx <= stop < d.end_tx and (stop - d.start_tx) % 3 == 0
                    for d in self.dorfs_by_transcript.get(transcript.transcript_id, [])
                ):
                    continue
                frame_start = stop % 3
                starts = [idx for idx in range(frame_start, stop, 3) if alt_seq[idx : idx + 3] in self.start_codons]
                starts = [idx for idx in starts if self._no_stop_between(alt_seq, idx + 3, stop)]
                if not starts:
                    continue
                for start in starts:
                    length = stop + 3 - start
                    if length < self.min_length or length > self.max_length:
                        continue
                    orf_class = self._classify_orf(transcript, start, stop + 3, alt_utr_start)
                    if not orf_class:
                        continue
                    bounds = self._alt_genomic_bounds(transcript, offset, delta, len(alt), start, stop + 3)
                    if bounds == (0, 0):
                        continue
                    start_codon = alt_seq[start : start + 3]
                    start_kozak = kozak_strength(alt_seq, start)
                    dummy = DORF(
                        dorf_id=f"{transcript.transcript_id}:{orf_class}:novel_stop:{variant.pos}:{variant.alt}:{start}",
                        transcript_id=transcript.transcript_id,
                        gene_id=transcript.gene_id,
                        gene_name=transcript.gene_name,
                        chrom=transcript.chrom,
                        strand=transcript.strand,
                        start_tx=start,
                        end_tx=stop + 3,
                        genomic_start=bounds[0],
                        genomic_end=bounds[1],
                        start_codon=start_codon,
                        stop_codon=alt_codon,
                        length_nt=length,
                        kozak_sequence=kozak_context(alt_seq, start),
                        kozak_strength=start_kozak,
                        orf_class=orf_class,
                        has_evidence=False,
                        transcript_dorf_count_total=self._evidence_orf_count(transcript.transcript_id),
                        dorf_start_distance_from_cds_stop_nt=start - alt_utr_start,
                        is_mane=transcript.is_mane,
                    )
                    effect = self._base_effect(dummy, variant)
                    effect.consequence = "dStop_gained"
                    effect.start_codon_alt = start_codon
                    effect.stop_codon_ref = ref_codon
                    effect.stop_codon_alt = alt_codon
                    effect.kozak_strength_alt = start_kozak
                    effect.has_evidence = False
                    effect.alt_dorf_length_nt = dummy.length_nt
                    stop_pos = self._alt_index_to_genomic(transcript, offset, delta, len(alt), stop)
                    effect.detail = f"new_stop={alt_codon};new_stop_pos={stop_pos};new_len={length};kozak={dummy.kozak_strength}"
                    effects.append(effect)
        return effects

    def _candidate_transcripts(self, variant: Variant) -> list[Transcript]:
        span_end = variant.pos + max(1, len(variant.ref)) - 1
        seen: set[str] = set()
        transcripts: list[Transcript] = []
        for alias in chrom_aliases(variant.chrom):
            for bin_id in self._bin_ids(variant):
                for transcript_id in self._transcript_bin_index.get((alias, bin_id), ()):
                    if transcript_id in seen:
                        continue
                    seen.add(transcript_id)
                    transcript = self.transcripts.get(transcript_id)
                    if transcript is None:
                        continue
                    if self.mane_only and not transcript.is_mane:
                        continue
                    bounds = self._transcript_bounds.get(transcript_id)
                    if bounds is not None and variant.pos <= bounds[1] and span_end >= bounds[0]:
                        transcripts.append(transcript)
        return transcripts

    def _evidence_orf_count(self, transcript_id: str) -> int:
        return sum(1 for dorf in self.dorfs_by_transcript.get(transcript_id, []) if dorf.has_evidence)

    def _next_stop(self, seq: str, start: int) -> int | None:
        for index in range(start, len(seq) - 2, 3):
            if seq[index : index + 3] in self.stop_codons:
                return index
        return None


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, coords = region.replace(",", "").split(":", 1)
    start, end = coords.split("-", 1)
    return chrom, int(start), int(end)


def open_variant_text(path: str | Path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def read_vcf(path: str | Path) -> Iterable[Variant]:
    skipped = 0
    with open_variant_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, pos, _vid, ref, alts, *_ = line.rstrip("\n").split("\t")
            for alt in alts.split(","):
                if alt == "." or not is_simple_allele(ref, alt):
                    skipped += 1
                    continue
                yield Variant(chrom, int(pos), ref, alt)
    if skipped:
        logger.warning("Skipped %d VCF allele(s) with symbolic, missing or non-ACGT sequence", skipped)


def read_tsv_variants(path: str | Path) -> Iterable[Variant]:
    skipped = 0
    with open_variant_text(path) as handle:
        field_map: dict[str, int] | None = None
        header_checked = False
        for line in handle:
            if not line.strip() or line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if not header_checked:
                normalized = [field.strip().lower().lstrip("#") for field in fields]
                if {"chrom", "pos", "ref", "alt"}.issubset(normalized):
                    field_map = {name: normalized.index(name) for name in ("chrom", "pos", "ref", "alt")}
                    header_checked = True
                    continue
                field_map = {"chrom": 0, "pos": 1, "ref": 2, "alt": 3}
                header_checked = True
            if line.startswith("#"):
                continue
            assert field_map is not None
            required_index = max(field_map.values())
            if len(fields) <= required_index:
                raise ValueError("Variant TSV rows must contain chrom, pos, ref, and alt columns")
            chrom = fields[field_map["chrom"]]
            pos = int(fields[field_map["pos"]])
            ref = fields[field_map["ref"]]
            alts = fields[field_map["alt"]]
            for alt in alts.split(","):
                alt = alt.strip()
                if not alt:
                    continue
                if not is_simple_allele(ref, alt):
                    skipped += 1
                    continue
                yield Variant(chrom, pos, ref, alt)
    if skipped:
        logger.warning("Skipped %d TSV allele(s) with non-ACGT sequence", skipped)



class EffectTsvWriter:
    def __init__(self, path: str | Path) -> None:
        self.path = path
        self.handle = open(path, "wt", encoding="utf-8", newline="")
        self.fieldnames = OUTPUT_FIELDS
        self.writer = csv.DictWriter(self.handle, fieldnames=self.fieldnames, delimiter="\t", lineterminator="\n")
        self.writer.writeheader()

    def write_many(self, effects: Iterable[Effect]) -> None:
        for effect in effects:
            self.writer.writerow(effect.to_row())

    def close(self) -> None:
        self.handle.close()

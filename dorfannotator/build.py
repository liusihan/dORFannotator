from __future__ import annotations

import gzip
import json
import logging
import sqlite3
from pathlib import Path
from typing import Iterable

logger = logging.getLogger(__name__)

from .models import DORF, Transcript
from .sequence import (
    DEFAULT_START_CODONS,
    STOP_CODONS,
    fetch_sequence_any,
    kozak_context,
    kozak_strength,
    normalize_chrom,
    read_fasta,
    revcomp,
)


def open_text(path: str | Path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def base_transcript_id(transcript_id: str) -> str:
    return transcript_id.split(".")[0]


def parse_gtf_attributes(raw: str) -> dict[str, str]:
    values: dict[str, str] = {}
    tags: list[str] = []
    for item in raw.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        clean = value.strip().strip('"')
        if key == "tag":
            tags.append(clean)
        else:
            values[key] = clean
    if tags:
        values["tag"] = ",".join(tags)
    return values


THREE_PRIME_UTR_FEATURES = {"three_prime_UTR", "three_prime_utr", "3UTR", "UTR3"}
GENERIC_UTR_FEATURES = {"UTR", "utr"}


def parse_gtf(path: str | Path) -> dict[str, dict]:
    transcripts: dict[str, dict] = {}
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs_raw = fields
            attrs = parse_gtf_attributes(attrs_raw)
            transcript_id = attrs.get("transcript_id")
            if feature == "transcript":
                transcript_id = transcript_id or attrs.get("ID")
            if not transcript_id:
                continue
            record = transcripts.setdefault(
                transcript_id,
                {
                    "transcript_id": transcript_id,
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_name": attrs.get("gene_name", attrs.get("gene_id", "")),
                    "chrom": chrom,
                    "strand": strand,
                    "utr3": [],
                    "utr_generic": [],
                    "cds": [],
                    "is_mane": False,
                },
            )
            record["gene_id"] = record["gene_id"] or attrs.get("gene_id", "")
            record["gene_name"] = record["gene_name"] or attrs.get("gene_name", record["gene_id"])
            tag_text = attrs.get("tag", "")
            if attrs.get("MANE_Select") or attrs.get("mane_select") or "MANE_Select" in tag_text or "MANE_Plus_Clinical" in tag_text:
                record["is_mane"] = True
            if feature in THREE_PRIME_UTR_FEATURES:
                record["utr3"].append((int(start), int(end)))
            elif feature in GENERIC_UTR_FEATURES:
                record["utr_generic"].append((int(start), int(end)))
            elif feature == "CDS":
                record["cds"].append((int(start), int(end)))
    for record in transcripts.values():
        generic = record.pop("utr_generic", [])
        if generic and record["cds"]:
            cds_max_end = max(end for _, end in record["cds"])
            cds_min_start = min(start for start, _ in record["cds"])
            for seg_start, seg_end in generic:
                if record["strand"] == "+":
                    if seg_start > cds_max_end:
                        record["utr3"].append((seg_start, seg_end))
                else:
                    if seg_end < cds_min_start:
                        record["utr3"].append((seg_start, seg_end))
        record["utr3"] = sorted(set(record["utr3"]))
    return transcripts


def build_transcript(record: dict, fasta: dict[str, str], max_length: int = 303) -> Transcript | None:
    if not record["utr3"]:
        return None
    chrom = record["chrom"]
    strand = record["strand"]
    segments = sorted(record["utr3"], key=lambda item: item[0], reverse=(strand == "-"))
    seq_parts: list[str] = []
    pos: list[int] = []
    for start, end in segments:
        segment = fetch_sequence_any(fasta, chrom, start, end)
        if strand == "-":
            seq_parts.append(revcomp(segment))
            pos.extend(range(end, start - 1, -1))
        else:
            seq_parts.append(segment)
            pos.extend(range(start, end + 1))
    cds_seq_parts: list[str] = []
    cds_pos: list[int] = []
    cds_segments = sorted(record["cds"], key=lambda item: item[0], reverse=(strand == "-"))
    for start, end in cds_segments:
        segment = fetch_sequence_any(fasta, chrom, start, end)
        if strand == "-":
            cds_seq_parts.append(revcomp(segment))
            cds_pos.extend(range(end, start - 1, -1))
        else:
            cds_seq_parts.append(segment)
            cds_pos.extend(range(start, end + 1))
    cds_seq = "".join(cds_seq_parts).upper()
    cds_tail_len = min(len(cds_seq), max_length)
    cds_tail = cds_seq[-cds_tail_len:] if cds_tail_len else ""
    cds_tail_pos = cds_pos[-cds_tail_len:] if cds_tail_len else []
    utr_sequence = "".join(seq_parts).upper()
    return Transcript(
        transcript_id=record["transcript_id"],
        gene_id=record["gene_id"],
        gene_name=record["gene_name"],
        chrom=chrom,
        strand=strand,
        utr_sequence=utr_sequence,
        utr_genomic_positions=tuple(pos),
        is_mane=bool(record.get("is_mane", False)),
        orf_sequence=cds_tail + utr_sequence,
        orf_genomic_positions=tuple(cds_tail_pos + pos),
        utr_start_tx=len(cds_tail),
        cds_tail_cds_start_offset=len(cds_seq) - cds_tail_len,
    )


def read_translated_dorf(path: str | Path | None) -> set[tuple[str, str, str, int, int, str]]:
    """Parse a translated dORF/doORF tab-delimited file.

    Required columns are:
        chrom  start  end  strand  transcript_id  orf_class
    """
    evidence: set[tuple[str, str, str, int, int, str]] = set()
    if not path:
        return evidence
    required = ("chrom", "start", "end", "strand", "transcript_id", "orf_class")
    field_map: dict[str, int] | None = None
    header_checked = False
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            normalized = [field.strip().lower().lstrip("#") for field in fields]
            if fields[0].startswith("#") and not set(required).issubset(normalized):
                continue
            if len(fields) < 6:
                raise ValueError(
                    "Translated dORF tab-delimited file must have at least 6 columns: "
                    "chrom  start  end  strand  transcript_id  orf_class"
                )
            if not header_checked:
                if set(required).issubset(normalized):
                    field_map = {name: normalized.index(name) for name in required}
                    header_checked = True
                    continue
                field_map = {name: index for index, name in enumerate(required)}
                header_checked = True
            assert field_map is not None
            try:
                start = int(fields[field_map["start"]])
                end = int(fields[field_map["end"]])
            except ValueError as exc:
                raise ValueError(
                    f"Translated dORF tab-delimited file line {line_number} start/end columns "
                    "must be integer 1-based genomic coordinates"
                ) from exc
            if start < 1 or end < 1:
                raise ValueError("Translated dORF tab-delimited file coordinates must be positive 1-based integers")
            if start > end:
                raise ValueError("Translated dORF tab-delimited file start must be less than or equal to end")
            strand = fields[field_map["strand"]]
            if strand not in {"+", "-"}:
                raise ValueError("Translated dORF tab-delimited file strand column must be '+' or '-'")
            transcript_id = base_transcript_id(fields[field_map["transcript_id"]])
            orf_class = fields[field_map["orf_class"]]
            if orf_class not in {"dORF", "doORF"}:
                raise ValueError("Translated dORF tab-delimited file orf_class column must be 'dORF' or 'doORF'")
            evidence.add((normalize_chrom(fields[field_map["chrom"]]), transcript_id, strand, start, end, orf_class))
    return evidence


def classify_orf(transcript: Transcript, start_tx: int, end_tx: int) -> str:
    stop_tx = end_tx - 3
    if start_tx >= transcript.utr_start_tx and stop_tx >= transcript.utr_start_tx:
        return "dORF"
    if start_tx < transcript.utr_start_tx <= stop_tx:
        if (transcript.cds_tail_cds_start_offset + start_tx) % 3 != 0:
            return "doORF"
    return ""


def scan_dorfs(
    transcript: Transcript,
    min_length: int,
    max_length: int,
    start_codons: Iterable[str],
    stop_codons: Iterable[str] = STOP_CODONS,
) -> list[DORF]:
    """Ab-initio scan for candidate dORFs/doORFs on a transcript.

    Every ORF returned here is sequence-predicted, so ``has_evidence`` is False;
    translated-evidence ORFs are produced separately by ``build_translated_dorfs``
    and combined in ``merge_dorfs``.
    """
    starts = {codon.upper() for codon in start_codons}
    stops = {codon.upper() for codon in stop_codons}
    seq = transcript.orf_sequence or transcript.utr_sequence
    positions = transcript.orf_genomic_positions or transcript.utr_genomic_positions
    dorfs: list[DORF] = []
    for start_tx in range(0, max(0, len(seq) - 2)):
        codon = seq[start_tx : start_tx + 3]
        if codon not in starts:
            continue
        for stop_tx in range(start_tx + 3, len(seq) - 2, 3):
            stop_codon = seq[stop_tx : stop_tx + 3]
            if stop_codon not in stops:
                continue
            length = stop_tx + 3 - start_tx
            orf_class = classify_orf(transcript, start_tx, stop_tx + 3)
            if not orf_class:
                break
            if min_length <= length <= max_length:
                genomic_positions = positions[start_tx : stop_tx + 3]
                genomic_start = min(genomic_positions)
                genomic_end = max(genomic_positions)
                strength = kozak_strength(seq, start_tx)
                dorf_id = f"{transcript.transcript_id}:{orf_class}:{genomic_start}-{genomic_end}"
                dorfs.append(
                    DORF(
                        dorf_id=dorf_id,
                        transcript_id=transcript.transcript_id,
                        gene_id=transcript.gene_id,
                        gene_name=transcript.gene_name,
                        chrom=transcript.chrom,
                        strand=transcript.strand,
                        start_tx=start_tx,
                        end_tx=stop_tx + 3,
                        genomic_start=genomic_start,
                        genomic_end=genomic_end,
                        start_codon=codon,
                        stop_codon=stop_codon,
                        length_nt=length,
                        kozak_sequence=kozak_context(seq, start_tx),
                        kozak_strength=strength,
                        orf_class=orf_class,
                        has_evidence=False,
                        dorf_start_distance_from_cds_stop_nt=start_tx - transcript.utr_start_tx,
                        is_mane=transcript.is_mane,
                    )
                )
            break
    dorfs.sort(key=lambda dorf: dorf.start_tx)
    total = len(dorfs)
    for dorf in dorfs:
        dorf.transcript_dorf_count_total = total
    return dorfs


def build_translated_dorfs(
    transcripts_by_tid: dict[str, list[Transcript]],
    translated_evidence: set[tuple[str, str, str, int, int, str]],
    stop_codons: Iterable[str] = STOP_CODONS,
) -> list[DORF]:
    """Build ORFs directly from translated dORF/doORF evidence.
    """
    stops = {codon.upper() for codon in stop_codons}
    index_cache: dict[str, dict[int, int]] = {}
    dorfs: list[DORF] = []
    for chrom, transcript_id, strand, gstart, gend, declared_class in translated_evidence:
        for transcript in transcripts_by_tid.get(transcript_id, []):
            if normalize_chrom(transcript.chrom) != chrom or transcript.strand != strand:
                continue
            pos_index = index_cache.get(transcript.transcript_id)
            if pos_index is None:
                positions = transcript.orf_genomic_positions or transcript.utr_genomic_positions
                pos_index = {pos: idx for idx, pos in enumerate(positions)}
                index_cache[transcript.transcript_id] = pos_index
            if gstart not in pos_index or gend not in pos_index:
                continue
            if strand == "+":
                start_tx, end_tx = pos_index[gstart], pos_index[gend] + 1
            else:
                start_tx, end_tx = pos_index[gend], pos_index[gstart] + 1
            seq = transcript.orf_sequence or transcript.utr_sequence
            if not (0 <= start_tx < end_tx <= len(seq)):
                continue
            sub = seq[start_tx:end_tx]
            if len(sub) < 3 or len(sub) % 3 != 0:
                continue
            if sub[-3:] not in stops:
                continue
            orf_class = declared_class
            strength = kozak_strength(seq, start_tx)
            dorfs.append(
                DORF(
                    dorf_id=f"{transcript.transcript_id}:{orf_class}:{gstart}-{gend}",
                    transcript_id=transcript.transcript_id,
                    gene_id=transcript.gene_id,
                    gene_name=transcript.gene_name,
                    chrom=transcript.chrom,
                    strand=transcript.strand,
                    start_tx=start_tx,
                    end_tx=end_tx,
                    genomic_start=gstart,
                    genomic_end=gend,
                    start_codon=sub[:3],
                    stop_codon=sub[-3:],
                    length_nt=len(sub),
                    kozak_sequence=kozak_context(seq, start_tx),
                    kozak_strength=strength,
                    orf_class=orf_class,
                    has_evidence=True,
                    dorf_start_distance_from_cds_stop_nt=start_tx - transcript.utr_start_tx,
                    is_mane=transcript.is_mane,
                )
            )
            break
    return dorfs


def _genomic_key(dorf: DORF) -> tuple[str, str, int, int, str]:
    return (base_transcript_id(dorf.transcript_id), dorf.strand,
            dorf.genomic_start, dorf.genomic_end, dorf.orf_class)


def merge_dorfs(translated: list[DORF], scanned: list[DORF]) -> list[DORF]:
    evidence_keys = {_genomic_key(d) for d in translated}
    extra = [d for d in scanned if _genomic_key(d) not in evidence_keys]
    return translated + extra


def finalize_transcript_dorf_stats(dorfs: list[DORF]) -> None:
    groups: dict[str, list[DORF]] = {}
    for dorf in dorfs:
        groups.setdefault(dorf.transcript_id, []).append(dorf)
    for group in groups.values():
        group.sort(key=lambda dorf: dorf.start_tx)
        total = sum(1 for dorf in group if dorf.has_evidence)
        for dorf in group:
            dorf.transcript_dorf_count_total = total


def init_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        DROP TABLE IF EXISTS transcripts;
        DROP TABLE IF EXISTS dorfs;
        DROP TABLE IF EXISTS dorf_bins;
        DROP TABLE IF EXISTS transcript_bins;
        DROP TABLE IF EXISTS metadata;
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        );
        CREATE TABLE transcripts (
            transcript_id TEXT PRIMARY KEY,
            gene_id TEXT,
            gene_name TEXT,
            chrom TEXT,
            strand TEXT,
            utr_sequence TEXT,
            utr_genomic_positions TEXT,
            is_mane INTEGER,
            orf_sequence TEXT,
            orf_genomic_positions TEXT,
            utr_start_tx INTEGER,
            cds_tail_cds_start_offset INTEGER
        );
        CREATE TABLE dorfs (
            dorf_id TEXT PRIMARY KEY,
            transcript_id TEXT,
            gene_id TEXT,
            gene_name TEXT,
            chrom TEXT,
            strand TEXT,
            start_tx INTEGER,
            end_tx INTEGER,
            genomic_start INTEGER,
            genomic_end INTEGER,
            start_codon TEXT,
            stop_codon TEXT,
            length_nt INTEGER,
            kozak_sequence TEXT,
            kozak_strength TEXT,
            orf_class TEXT,
            has_evidence INTEGER,
            transcript_dorf_count_total INTEGER,
            dorf_start_distance_from_cds_stop_nt INTEGER,
            is_mane INTEGER
        );
        CREATE TABLE dorf_bins (
            chrom TEXT,
            bin INTEGER,
            dorf_id TEXT
        );
        CREATE TABLE transcript_bins (
            chrom TEXT,
            bin INTEGER,
            transcript_id TEXT
        );
        CREATE INDEX idx_dorf_bins ON dorf_bins(chrom, bin);
        CREATE INDEX idx_transcript_bins ON transcript_bins(chrom, bin);
        CREATE INDEX idx_dorfs_transcript ON dorfs(transcript_id);
        """
    )


def write_database(
    db_path: str | Path,
    transcripts: list[Transcript],
    dorfs: list[DORF],
    bin_size: int,
    min_length: int,
    max_length: int,
    start_codons: Iterable[str],
    stop_codons: Iterable[str],
) -> None:
    ensure_unique_dorf_ids(dorfs)
    conn = sqlite3.connect(db_path)
    try:
        init_schema(conn)
        conn.executemany(
            "INSERT INTO metadata VALUES (?, ?)",
            [
                ("min_length", str(min_length)),
                ("max_length", str(max_length)),
                ("bin_size", str(bin_size)),
                ("start_codons", ",".join(start_codons)),
                ("stop_codons", ",".join(stop_codons)),
            ],
        )
        conn.executemany(
            """
            INSERT INTO transcripts VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    tx.transcript_id,
                    tx.gene_id,
                    tx.gene_name,
                    tx.chrom,
                    tx.strand,
                    tx.utr_sequence,
                    json.dumps(tx.utr_genomic_positions),
                    1 if tx.is_mane else 0,
                    tx.orf_sequence,
                    json.dumps(tx.orf_genomic_positions),
                    tx.utr_start_tx,
                    tx.cds_tail_cds_start_offset,
                )
                for tx in transcripts
            ],
        )
        conn.executemany(
            """
            INSERT INTO dorfs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    d.dorf_id,
                    d.transcript_id,
                    d.gene_id,
                    d.gene_name,
                    d.chrom,
                    d.strand,
                    d.start_tx,
                    d.end_tx,
                    d.genomic_start,
                    d.genomic_end,
                    d.start_codon,
                    d.stop_codon,
                    d.length_nt,
                    d.kozak_sequence,
                    d.kozak_strength,
                    d.orf_class,
                    1 if d.has_evidence else 0,
                    d.transcript_dorf_count_total,
                    d.dorf_start_distance_from_cds_stop_nt,
                    1 if d.is_mane else 0,
                )
                for d in dorfs
            ],
        )
        bin_rows = []
        for dorf in dorfs:
            pad_low = 6 if dorf.strand == "+" else 0
            pad_high = 6 if dorf.strand == "-" else 0
            padded_start = max(1, dorf.genomic_start - pad_low)
            padded_end = dorf.genomic_end + pad_high
            for bin_id in range(padded_start // bin_size, padded_end // bin_size + 1):
                bin_rows.append((dorf.chrom, bin_id, dorf.dorf_id))
        conn.executemany("INSERT INTO dorf_bins VALUES (?, ?, ?)", bin_rows)
        transcript_bin_rows = []
        for tx in transcripts:
            positions = tx.orf_genomic_positions or tx.utr_genomic_positions
            if not positions:
                continue
            start = min(positions)
            end = max(positions)
            for bin_id in range(start // bin_size, end // bin_size + 1):
                transcript_bin_rows.append((tx.chrom, bin_id, tx.transcript_id))
        conn.executemany("INSERT INTO transcript_bins VALUES (?, ?, ?)", transcript_bin_rows)
        conn.commit()
    finally:
        conn.close()


def ensure_unique_dorf_ids(dorfs: list[DORF]) -> None:
    seen: dict[str, int] = {}
    for dorf in dorfs:
        base_id = dorf.dorf_id
        seen[base_id] = seen.get(base_id, 0) + 1
        if seen[base_id] > 1:
            dorf.dorf_id = f"{base_id}#{seen[base_id]}"


def read_transcript_id_file(path: str | Path) -> set[str]:
    ids: set[str] = set()
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.add(line.split()[0])
    return ids


def build_database(
    gtf_path: str | Path,
    fasta_path: str | Path,
    translated_dorf_path: str | Path | None,
    db_path: str | Path,
    min_length: int = 30,
    max_length: int = 303,
    start_codons: Iterable[str] = DEFAULT_START_CODONS,
    stop_codons: Iterable[str] = STOP_CODONS,
    bin_size: int = 10000,
    mane_only: bool = False,
    mane_transcripts: str | Path | None = None,
) -> None:
    logger.info("Loading reference FASTA %s", fasta_path)
    fasta = read_fasta(str(fasta_path))
    logger.info("Parsing GTF %s", gtf_path)
    records = parse_gtf(gtf_path)
    logger.info("Parsed %d transcripts from GTF", len(records))
    translated = read_translated_dorf(translated_dorf_path)
    logger.info("Read %d translated-evidence entries", len(translated))
    mane_ids = {base_transcript_id(tid) for tid in read_transcript_id_file(mane_transcripts)} if mane_transcripts else set()

    def apply_mane(tx: Transcript) -> Transcript:
        if mane_ids and base_transcript_id(tx.transcript_id) in mane_ids:
            return Transcript(**{**tx.__dict__, "is_mane": True})
        return tx

    FULL_CDS_CONTEXT = 10 ** 9  # include the entire CDS so deep doORF starts stay representable
    evidence_tids = {tid for (_c, tid, _s, _gs, _ge, _cls) in translated}
    evidence_by_tid: dict[str, list[Transcript]] = {}
    evidence_tx_by_id: dict[str, Transcript] = {}
    for record in records.values():
        if base_transcript_id(record["transcript_id"]) not in evidence_tids:
            continue
        tx = build_transcript(record, fasta, max_length=FULL_CDS_CONTEXT)
        if tx is None:
            continue
        tx = apply_mane(tx)
        if mane_only and not tx.is_mane:
            continue
        evidence_by_tid.setdefault(base_transcript_id(tx.transcript_id), []).append(tx)
        evidence_tx_by_id[tx.transcript_id] = tx
    translated_dorfs = build_translated_dorfs(evidence_by_tid, translated, stop_codons)
    evidence_keys = {_genomic_key(d) for d in translated_dorfs}
    logger.info(
        "Phase 1: matched %d translated-evidence ORFs from %d input entries",
        len(translated_dorfs), len(translated),
    )

    transcripts: list[Transcript] = []
    predicted: list[DORF] = []
    for record in records.values():
        full_tx = evidence_tx_by_id.get(record["transcript_id"])
        if full_tx is not None:
            tx = full_tx
        else:
            tx = build_transcript(record, fasta, max_length=max_length)
            if tx is None:
                continue
            tx = apply_mane(tx)
            if mane_only and not tx.is_mane:
                continue
        transcripts.append(tx)
        for dorf in scan_dorfs(tx, min_length, max_length, start_codons, stop_codons):
            if _genomic_key(dorf) not in evidence_keys:
                predicted.append(dorf)

    logger.info(
        "Phase 2: predicted %d ab-initio ORFs over %d transcripts",
        len(predicted), len(transcripts),
    )
    dorfs = merge_dorfs(translated_dorfs, predicted)
    finalize_transcript_dorf_stats(dorfs)
    logger.info("Writing database %s (%d transcripts, %d ORFs)", db_path, len(transcripts), len(dorfs))
    write_database(db_path, transcripts, dorfs, bin_size, min_length, max_length, tuple(start_codons), tuple(stop_codons))
    logger.info("Build complete")

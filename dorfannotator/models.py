from __future__ import annotations

from dataclasses import dataclass
from typing import Any


OUTPUT_FIELDS = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "gene",
    "transcript",
    "mane",
    "strand",
    "orf_class",
    "csq",
    "evidence",
    "dorf_count",
    "dorf_start",
    "dorf_end",
    "dist_cds",
    "detail",
]


@dataclass(frozen=True)
class Transcript:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    utr_sequence: str
    utr_genomic_positions: tuple[int, ...]
    is_mane: bool = False
    orf_sequence: str = ""
    orf_genomic_positions: tuple[int, ...] = ()
    utr_start_tx: int = 0
    cds_tail_cds_start_offset: int = 0


@dataclass
class DORF:
    dorf_id: str
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    start_tx: int
    end_tx: int
    genomic_start: int
    genomic_end: int
    start_codon: str
    stop_codon: str
    length_nt: int
    kozak_sequence: str
    kozak_strength: str
    orf_class: str = "dORF"
    has_evidence: bool = False
    transcript_dorf_count_total: int = 0
    dorf_start_distance_from_cds_stop_nt: int = 0
    is_mane: bool = False


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass
class Effect:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene_id: str
    gene_name: str
    transcript_id: str
    is_mane: bool
    strand: str
    dorf_id: str
    orf_class: str
    consequence: str
    transcript_dorf_count_total: int = 0
    dorf_start_distance_from_cds_stop_nt: int = 0
    start_codon_ref: str = ""
    start_codon_alt: str = ""
    stop_codon_ref: str = ""
    stop_codon_alt: str = ""
    kozak_strength_ref: str = ""
    kozak_strength_alt: str = ""
    ref_dorf_length_nt: int = 0
    alt_dorf_length_nt: int = 0
    dorf_length_delta_nt: int = 0
    alternative_start_found: bool = False
    alternative_start_position: int | None = None
    alternative_start_kozak_strength: str = ""
    alternative_start_codon: str = ""
    alternative_stop_found: bool = False
    alternative_stop_position: int | None = None
    alternative_stop_codon: str = ""
    extension_length_nt: int = 0
    truncation_length_nt: int = 0
    has_evidence: bool = False
    dorf_start: int = 0
    dorf_end: int = 0
    dist_cds: int = 0
    detail: str = ""

    def to_row(self) -> dict[str, Any]:
        row: dict[str, Any] = {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "gene": self.gene_name or self.gene_id,
            "transcript": self.transcript_id,
            "mane": self.is_mane,
            "strand": self.strand,
            "orf_class": self.orf_class,
            "csq": self.consequence,
            "evidence": self.has_evidence,
            "dorf_count": self.transcript_dorf_count_total,
            "dorf_start": self.dorf_start,
            "dorf_end": self.dorf_end,
            "dist_cds": self.dist_cds,
            "detail": self.detail,
        }
        for key, value in list(row.items()):
            if isinstance(value, bool):
                row[key] = "true" if value else "false"
            elif value is None:
                row[key] = ""
            else:
                row[key] = str(value)
        return row

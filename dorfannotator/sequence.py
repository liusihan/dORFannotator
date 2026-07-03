from __future__ import annotations

import gzip


DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOP_CODONS = {"TAA", "TAG", "TGA"}
DEFAULT_START_CODONS = ("ATG", "CTG", "TTG", "GTG", "AAG", "ACG", "AGG", "ATC", "ATA", "ATT")


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMPLEMENT)[::-1].upper()


def read_fasta(path: str) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                sequences[current] = []
            elif current is not None:
                sequences[current].append(line.upper())
    return {chrom: "".join(parts) for chrom, parts in sequences.items()}


def fetch_sequence(fasta: dict[str, str], chrom: str, start: int, end: int) -> str:
    return fasta[chrom][start - 1 : end].upper()


def kozak_context(seq: str, start_index: int) -> str:
    left = max(0, start_index - 6)
    right = min(len(seq), start_index + 7)
    return seq[left:right]


def kozak_strength(seq: str, start_index: int) -> str:
    minus3 = seq[start_index - 3].upper() if start_index >= 3 else "N"
    plus4_index = start_index + 3
    plus4 = seq[plus4_index].upper() if plus4_index < len(seq) else "N"
    minus3_ok = minus3 in {"A", "G"}
    plus4_ok = plus4 == "G"
    if minus3_ok and plus4_ok:
        return "strong"
    if minus3_ok or plus4_ok:
        return "moderate"
    return "weak"


def apply_variant_to_sequence(seq: str, offset: int, ref: str, alt: str) -> str:
    return seq[:offset] + alt.upper() + seq[offset + len(ref) :]



def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def chrom_aliases(chrom: str) -> tuple[str, ...]:
    norm = normalize_chrom(chrom)
    aliases = [chrom]
    if norm not in aliases:
        aliases.append(norm)
    prefixed = f"chr{norm}"
    if prefixed not in aliases:
        aliases.append(prefixed)
    return tuple(aliases)


def fetch_sequence_any(fasta: dict[str, str], chrom: str, start: int, end: int) -> str:
    for alias in chrom_aliases(chrom):
        if alias in fasta:
            return fetch_sequence(fasta, alias, start, end)
    raise KeyError(f"Chromosome {chrom!r} was not found in FASTA. Tried aliases: {', '.join(chrom_aliases(chrom))}")

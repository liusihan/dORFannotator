from __future__ import annotations

import argparse
import logging
from concurrent.futures import ProcessPoolExecutor

from . import __version__
from .annotator import Annotator, EffectTsvWriter, read_tsv_variants, read_vcf
from .build import build_database


def annotate_worker(db_path: str, variants, chromosomes, region, mane_only, include_predicted, evidence_only):
    load_chromosomes = sorted({variant.chrom for variant in variants})
    annotator = Annotator(
        db_path,
        chromosomes=chromosomes,
        region=region,
        load_chromosomes=load_chromosomes,
        mane_only=mane_only,
        include_predicted=include_predicted,
        evidence_only=evidence_only,
    )
    try:
        return annotator.annotate_batch(variants)
    finally:
        annotator.close()


def split_csv(value: str | None) -> tuple[str, ...]:
    if not value:
        return tuple()
    return tuple(item.strip() for item in value.split(",") if item.strip())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="dORFannotator")
    parser.add_argument("--version", action="version", version=f"dORFannotator {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build", help="Build a dORF/doORF annotation database")
    build.add_argument("--gtf", required=True)
    build.add_argument("--fasta", required=True)
    build.add_argument(
        "--translated-dorf",
        default=None,
        help="Optional tab-delimited translated dORF/doORF file with columns: chrom, start, end, strand, transcript_id, orf_class",
    )
    build.add_argument("--out-db", required=True)
    build.add_argument("--start-codons", default="ATG,CTG,TTG,GTG,AAG,ACG,AGG,ATC,ATA,ATT")
    build.add_argument("--stop-codons", default="TAG,TAA,TGA")
    build.add_argument("--min-length", type=int, default=30)
    build.add_argument("--max-length", type=int, default=303)
    build.add_argument("--mane-only", action="store_true", help="Build only transcripts tagged as MANE in the GTF")
    build.add_argument("--mane-transcripts", help="File containing MANE transcript IDs, one per line")
    build.add_argument("--verbose", action="store_true")

    annotate = subparsers.add_parser("annotate", help="Annotate VCF or TSV variants")
    variant_input = annotate.add_mutually_exclusive_group(required=True)
    variant_input.add_argument("--vcf", help="VCF/VCF.GZ variant input")
    variant_input.add_argument("--tsv", help="TSV/TSV.GZ variant input with chrom, pos, ref, alt columns")
    annotate.add_argument("--db", required=True)
    annotate.add_argument("--out", required=True)
    annotate.add_argument("--chr", dest="chromosomes", help="Comma-separated chromosome filter, e.g. chr1,chr2,chrX")
    annotate.add_argument("--region", help="Genomic region filter, e.g. chr1:1000000-2000000")
    annotate.add_argument("--threads", type=int, default=1)
    annotate.add_argument("--batch-size", type=int, default=10000)
    annotate.add_argument("--mane-only", action="store_true", help="Only output MANE transcript annotations")
    output_mode = annotate.add_mutually_exclusive_group()
    output_mode.add_argument("--include-predicted", action="store_true", help="Include all sequence-predicted ORF consequences")
    output_mode.add_argument("--evidence-only", action="store_true", help="Only output consequences on translated-evidence ORFs")
    annotate.add_argument("--verbose", action="store_true")

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO if getattr(args, "verbose", False) else logging.WARNING, format="%(levelname)s: %(message)s")

    if args.command == "build":
        build_database(
            gtf_path=args.gtf,
            fasta_path=args.fasta,
            translated_dorf_path=args.translated_dorf,
            db_path=args.out_db,
            min_length=args.min_length,
            max_length=args.max_length,
            start_codons=split_csv(args.start_codons),
            stop_codons=split_csv(args.stop_codons),
            mane_only=args.mane_only,
            mane_transcripts=args.mane_transcripts,
        )
        return 0

    chromosomes = split_csv(args.chromosomes)
    variant_reader = read_vcf if args.vcf else read_tsv_variants
    variant_path = args.vcf or args.tsv
    writer = EffectTsvWriter(args.out)
    try:
        if args.threads <= 1:
            annotator = Annotator(
                args.db,
                chromosomes=chromosomes,
                region=args.region,
                mane_only=args.mane_only,
                include_predicted=args.include_predicted,
                evidence_only=args.evidence_only,
            )
            try:
                batch = []
                for variant in variant_reader(variant_path):
                    batch.append(variant)
                    if len(batch) >= args.batch_size:
                        writer.write_many(annotator.annotate_batch(batch))
                        batch.clear()
                if batch:
                    writer.write_many(annotator.annotate_batch(batch))
            finally:
                annotator.close()
        else:
            pending = []
            with ProcessPoolExecutor(max_workers=args.threads) as pool:
                batch = []
                for variant in variant_reader(variant_path):
                    batch.append(variant)
                    if len(batch) >= args.batch_size:
                        pending.append(
                            pool.submit(
                                annotate_worker,
                                args.db,
                                batch,
                                chromosomes,
                                args.region,
                                args.mane_only,
                                args.include_predicted,
                                args.evidence_only,
                            )
                        )
                        batch = []
                    while len(pending) >= args.threads * 2:
                        writer.write_many(pending.pop(0).result())
                if batch:
                    pending.append(
                        pool.submit(
                            annotate_worker,
                            args.db,
                            batch,
                            chromosomes,
                            args.region,
                            args.mane_only,
                            args.include_predicted,
                            args.evidence_only,
                        )
                    )
                for future in pending:
                    writer.write_many(future.result())
    finally:
        writer.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

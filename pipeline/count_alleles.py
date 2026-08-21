#!/usr/bin/env python3
"""Count REF/ALT bases at donor heterozygous SNVs in one filtered BAM."""
import argparse
import csv
import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-baseq", type=int, default=0)
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    vcf = pysam.VariantFile(args.vcf)
    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["chrom", "pos", "id", "ref", "alt", "ref_count", "alt_count", "other_count", "depth"])
        for rec in vcf.fetch():
            if len(rec.ref) != 1 or len(rec.alts or ()) != 1 or len(rec.alts[0]) != 1:
                continue
            counts = {base: 0 for base in "ACGTN"}
            for col in bam.pileup(rec.contig, rec.pos - 1, rec.pos, truncate=True,
                                  stepper="samtools", min_base_quality=args.min_baseq,
                                  max_depth=10000):
                if col.reference_pos != rec.pos - 1:
                    continue
                for pileup_read in col.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue
                    base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                    counts[base if base in counts else "N"] += 1
            ref_n = counts.get(rec.ref.upper(), 0)
            alt_n = counts.get(rec.alts[0].upper(), 0)
            depth = sum(counts.values())
            writer.writerow([rec.contig, rec.pos, rec.id or ".", rec.ref, rec.alts[0],
                             ref_n, alt_n, depth - ref_n - alt_n, depth])


if __name__ == "__main__":
    main()


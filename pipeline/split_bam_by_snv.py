#!/usr/bin/env python3
"""Very simple SNV-based BAM splitter for a training-data prototype.

The two output BAMs retain the original reference coordinates. Reads are
assigned only when their observed bases support one of the two alleles at an
SNV; reads with no usable or conflicting evidence are written to ambiguous.tsv.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--snvs", required=True, help="TSV with chrom and pos; optionally ref and alt")
    ap.add_argument("--ref-bam", required=True)
    ap.add_argument("--alt-bam", required=True)
    ap.add_argument("--ambiguous", required=True)
    ap.add_argument("--window", type=int, default=500)
    args = ap.parse_args()

    snvs = defaultdict(list)
    with open(args.snvs, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom = row.get("chrom") or row.get("#chr")
            pos = int(row.get("pos") or row.get("end"))
            ref = (row.get("ref") or "").upper()
            alt = (row.get("alt") or "").upper()
            if chrom and ref and alt:
                snvs[chrom].append((pos, ref, alt))
    for chrom in snvs:
        snvs[chrom].sort()

    bam = pysam.AlignmentFile(args.bam, "rb")
    ref_bam = pysam.AlignmentFile(args.ref_bam, "wb", template=bam)
    alt_bam = pysam.AlignmentFile(args.alt_bam, "wb", template=bam)
    assigned = {}
    with open(args.ambiguous, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["read_name", "status", "evidence"])
        for chrom, records in snvs.items():
            for pos, ref, alt in records:
                for read in bam.fetch(chrom, max(0, pos - 1 - args.window), pos - 1 + args.window):
                    if read.is_unmapped or read.query_sequence is None:
                        continue
                    observed = None
                    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                        if ref_pos == pos - 1:
                            observed = read.query_sequence[query_pos].upper()
                            break
                    if observed not in (ref, alt):
                        continue
                    target = "ref" if observed == ref else "alt"
                    name = read.query_name
                    if name in assigned and assigned[name] != target:
                        assigned[name] = "ambiguous"
                        writer.writerow([name, "conflicting", f"{assigned[name]}->{target}"])
                    else:
                        assigned[name] = target

    bam.reset()
    for read in bam.fetch(until_eof=True):
        target = assigned.get(read.query_name)
        if target == "ref":
            ref_bam.write(read)
        elif target == "alt":
            alt_bam.write(read)
        elif target == "ambiguous":
            continue
    bam.close(); ref_bam.close(); alt_bam.close()
    pysam.index(args.ref_bam)
    pysam.index(args.alt_bam)


if __name__ == "__main__":
    main()

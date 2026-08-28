#!/usr/bin/env python3
"""Split one chromosome's reads into two synthetic allele BAMs."""
import argparse
import random
from collections import defaultdict
from pathlib import Path

import pysam
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--snvs", required=True)
    ap.add_argument("--chrom", default="chr1")
    ap.add_argument("--output-prefix", required=True)
    ap.add_argument("--window", type=int, default=300)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    table = pd.read_csv(args.snvs, sep="\t")
    table = table[table["chrom"].astype(str).isin([args.chrom, args.chrom.removeprefix("chr")])]
    snvs = {int(r.pos): (r.ref.upper(), r.alt.upper()) for r in table.itertuples()}
    allele = {p: ((ref, alt) if rng.random() < .5 else (alt, ref))
              for p, (ref, alt) in snvs.items()}

    bam = pysam.AlignmentFile(args.bam, "rb")
    assignments = {}
    evidence = defaultdict(list)

    # Assign reads carrying an informative SNV.
    for pos, (a, b) in allele.items():
        for read in bam.fetch(args.chrom, pos - 1, pos):
            if read.is_unmapped or read.query_sequence is None:
                continue
            observed = None
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos == pos - 1:
                    observed = read.query_sequence[qpos].upper()
                    break
            if observed not in (a, b):
                continue
            key = (read.query_name, read.is_read1, read.is_read2)
            target = "A" if observed == a else "B"
            evidence[key].append((pos, target))
            if key not in assignments:
                assignments[key] = target
            elif assignments[key] != target:
                assignments[key] = "AMBIGUOUS"

    # Estimate local proportions from informative reads.
    local = defaultdict(lambda: {"A": 0, "B": 0})
    for key, obs in evidence.items():
        if assignments.get(key) in ("A", "B"):
            for pos, target in obs:
                local[pos][target] += 1

    # Assign every remaining mapped read: local proportion near an SNV, else 50/50.
    for read in bam.fetch(args.chrom):
        if read.is_unmapped:
            continue
        key = (read.query_name, read.is_read1, read.is_read2)
        if key in assignments:
            continue
        nearby = [p for p in allele if read.reference_start <= p - 1 + args.window
                   and read.reference_end >= p - args.window]
        a = sum(local[p]["A"] for p in nearby)
        b = sum(local[p]["B"] for p in nearby)
        probability_a = a / (a + b) if a + b else 0.5
        assignments[key] = "A" if rng.random() < probability_a else "B"

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    out_a = str(prefix) + ".A.bam"
    out_b = str(prefix) + ".B.bam"
    bam.reset()
    a_bam = pysam.AlignmentFile(out_a, "wb", template=bam)
    b_bam = pysam.AlignmentFile(out_b, "wb", template=bam)
    for read in bam.fetch(args.chrom):
        key = (read.query_name, read.is_read1, read.is_read2)
        if assignments.get(key) == "A":
            a_bam.write(read)
        elif assignments.get(key) == "B":
            b_bam.write(read)
    a_bam.close(); b_bam.close(); bam.close()
    pysam.index(out_a); pysam.index(out_b)
    print(f"Wrote {out_a} and {out_b}")


if __name__ == "__main__":
    main()

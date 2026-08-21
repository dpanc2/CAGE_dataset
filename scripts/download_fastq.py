#!/usr/bin/env python3
"""Download selected FASTQs from a manifest and verify their ENA MD5 values."""

import argparse
import csv
import hashlib
import subprocess
from pathlib import Path


def md5sum(path: Path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--runs", nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("fastq"))
    args = parser.parse_args()

    with args.manifest.open(newline="") as handle:
        rows = {r["run_accession"]: r for r in csv.DictReader(handle, delimiter="\t")}
    missing = sorted(set(args.runs) - rows.keys())
    if missing:
        raise SystemExit(f"Runs absent from manifest: {', '.join(missing)}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for accession in args.runs:
        row = rows[accession]
        destination = args.output_dir / f"{accession}.fastq.gz"
        partial = destination.with_suffix(destination.suffix + ".part")
        print(f"Downloading {accession} ({row['fastq_gib']} GiB)...", flush=True)
        subprocess.run(
            ["curl", "-L", "--fail", "--show-error", "--continue-at", "-",
             "--retry", "10", "--retry-all-errors", "--retry-delay", "2",
             "--speed-limit", "1024", "--speed-time", "30",
             "--output", str(partial), row["fastq_url"]],
            check=True,
        )
        observed = md5sum(partial)
        if observed != row["fastq_md5"]:
            raise SystemExit(f"MD5 mismatch for {accession}: {observed}")
        partial.replace(destination)
        print(f"Verified {destination}", flush=True)


if __name__ == "__main__":
    main()

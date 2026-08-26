#!/usr/bin/env python3
"""Restartable ENA FASTQ downloader driven by the all-sample manifest."""

import argparse
import csv
import hashlib
import subprocess
from pathlib import Path


def md5sum(path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with args.samples.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    for number, row in enumerate(rows, start=1):
        run = row["run_accession"]
        final = args.output_dir / row["fastq_name"]
        partial = Path(f"{final}.part")
        expected_size = int(row["fastq_bytes"])
        expected_md5 = row["fastq_md5"]

        if final.is_file() and final.stat().st_size == expected_size:
            print(f"[{number}/{len(rows)}] already present: {run}", flush=True)
            continue
        if final.exists():
            raise SystemExit(f"Existing file has wrong size: {final}")

        print(f"[{number}/{len(rows)}] downloading {run}", flush=True)
        subprocess.run([
            "curl", "-L", "--fail", "--show-error", "--continue-at", "-",
            "--retry", "10", "--retry-all-errors", "--retry-delay", "2",
            "--speed-limit", "1024", "--speed-time", "30",
            "--output", str(partial), row["fastq_url"],
        ], check=True)
        if partial.stat().st_size != expected_size:
            raise SystemExit(f"Size mismatch for {run}")
        observed_md5 = md5sum(partial)
        if observed_md5 != expected_md5:
            raise SystemExit(f"MD5 mismatch for {run}: {observed_md5}")
        partial.replace(final)
        print(f"[{number}/{len(rows)}] verified {run}", flush=True)


if __name__ == "__main__":
    main()

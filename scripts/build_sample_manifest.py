#!/usr/bin/env python3
"""Join GEO sample annotations to ENA run/FASTQ metadata for GSE150736."""

import argparse
import csv
import gzip
import re
from collections import Counter
from pathlib import Path


def parse_geo_soft(path: Path):
    samples = []
    current = None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    samples.append(current)
                current = {"geo_accession": line.split(" = ", 1)[1]}
            elif current is not None and line.startswith("!Sample_title = "):
                current["sample_title"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                key, value = value.split(":", 1)
                current[key.strip()] = value.strip()
            elif current is not None and line.startswith("!Sample_relation = BioSample: "):
                current["sample_accession"] = line.rsplit("/", 1)[-1]
            elif current is not None and line.startswith("!Sample_relation = SRA: "):
                match = re.search(r"term=(SRX\d+)", line)
                if match:
                    current["experiment_accession"] = match.group(1)
        if current:
            samples.append(current)
    return samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--soft", type=Path, required=True)
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    geo = parse_geo_soft(args.soft)
    with args.runs.open(newline="") as handle:
        runs = list(csv.DictReader(handle, delimiter="\t"))
    run_by_experiment = {row["experiment_accession"]: row for row in runs}

    tuples = [(x.get("gender", ""), x.get("age", ""), x.get("state", "")) for x in geo]
    tuple_counts = Counter(tuples)
    fields = [
        "geo_accession", "experiment_accession", "run_accession", "sample_accession",
        "sample_title", "protocol", "gender", "age", "state", "chamber",
        "donor_key_inferred", "donor_key_count", "donor_pair_confidence",
        "library_layout", "fastq_url", "fastq_md5", "fastq_bytes", "fastq_gib",
    ]
    output = []
    for sample, donor_tuple in zip(geo, tuples):
        run = run_by_experiment.get(sample.get("experiment_accession", ""), {})
        count = tuple_counts[donor_tuple]
        protocol = sample.get("sample_title", "").split(",", 1)[0]
        confidence = "candidate_unique_pair" if count == 2 else "ambiguous"
        size = int((run.get("fastq_bytes") or "0").split(";")[0])
        ftp = run.get("fastq_ftp", "").split(";")[0]
        output.append({
            **{key: sample.get(key, "") for key in fields},
            "run_accession": run.get("run_accession", ""),
            "sample_accession": run.get("sample_accession", sample.get("sample_accession", "")),
            "protocol": protocol,
            "donor_key_inferred": "|".join(donor_tuple),
            "donor_key_count": count,
            "donor_pair_confidence": confidence,
            "library_layout": run.get("library_layout", ""),
            "fastq_url": f"https://{ftp}" if ftp else "",
            "fastq_md5": (run.get("fastq_md5") or "").split(";")[0],
            "fastq_bytes": size,
            "fastq_gib": f"{size / 1024**3:.3f}",
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()

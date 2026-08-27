#!/usr/bin/env bash
set -euo pipefail

# Restartable preprocessing of all 109 GSE150736 CAGE libraries.
# STAGE may be: prepare, download, align, call, count, or all.

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${CAGE_DATA_ROOT:-/mnt/newdata/dpanc/CAGE_dataset_data}"
THREADS="${THREADS:-8}"
STAGE="${STAGE:-all}"
GEO_MANIFEST="${GEO_MANIFEST:-$REPO_DIR/data/metadata/GSE150736_sample_manifest.tsv}"
META_DIR="$DATA_ROOT/reference/heart_atlas_metadata"
SUPPLEMENT="$META_DIR/Supplementary_Table_1.xlsx"
SAMPLES="$META_DIR/GSE150736_all_samples.tsv"
FASTQ_DIR="$DATA_ROOT/fastq"
REF_DIR="$DATA_ROOT/reference/gencode_v39"
RESULTS="$DATA_ROOT/results/all"
FA="$REF_DIR/GRCh38.primary_assembly.genome.fa"
INDEX="$REF_DIR/hisat2/GRCh38_gencode_v39"

need() { command -v "$1" >/dev/null || { echo "Missing required command: $1" >&2; exit 1; }; }
for cmd in curl hisat2 samtools bcftools python3; do need "$cmd"; done

mkdir -p "$META_DIR" "$FASTQ_DIR" "$RESULTS"/{bam,logs,variants,counts,tmp}

run_stage() {
  [[ "$STAGE" == all || "$STAGE" == "$1" ]]
}

prepare() {
  [[ -s "$GEO_MANIFEST" ]] || { echo "Missing manifest: $GEO_MANIFEST" >&2; exit 1; }
  if [[ ! -s "$SUPPLEMENT" ]]; then
    curl -L --fail --retry 10 \
      'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs44161-022-00182-x/MediaObjects/44161_2022_182_MOESM3_ESM.xlsx' \
      -o "$SUPPLEMENT.part"
    mv "$SUPPLEMENT.part" "$SUPPLEMENT"
  fi
  python3 "$REPO_DIR/pipeline/build_all_samples.py" \
    --manifest "$GEO_MANIFEST" --supplement "$SUPPLEMENT" --output "$SAMPLES"
}

if run_stage prepare || [[ ! -s "$SAMPLES" ]]; then prepare; fi

if run_stage download; then
  python3 "$REPO_DIR/pipeline/download_all_fastq.py" \
    --samples "$SAMPLES" --output-dir "$FASTQ_DIR"
fi

if run_stage align; then
  awk -F '\t' 'NR>1 {print $1 "\t" $3 "\t" $6 "\t" $11}' "$SAMPLES" \
  | while IFS=$'\t' read -r donor run chamber fq_name; do
    fq="$FASTQ_DIR/$fq_name"
    out="$RESULTS/bam/$run.filtered.bam"
    [[ -s "$fq" ]] || { echo "Missing FASTQ: $fq" >&2; exit 1; }
    if [[ -s "$out" && -s "$out.bai" ]]; then
      echo "Alignment already complete: $run"
      continue
    fi
    echo "Aligning $run ($donor, $chamber)"
    tmp_prefix="$RESULTS/tmp/$run"
    hisat2 -p "$THREADS" --very-sensitive -x "$INDEX" -U "$fq" \
      --rg-id "$run" --rg "SM:$donor" \
      --summary-file "$RESULTS/logs/$run.hisat2.txt" \
      2> "$RESULTS/logs/$run.stderr.txt" \
    | samtools view -@ "$THREADS" -u - \
    | samtools sort -@ "$THREADS" -m 1G -T "$tmp_prefix.namesort" -n \
        -o "$RESULTS/tmp/$run.namesort.bam" -

    python3 "$REPO_DIR/pipeline/vendor/filter_reads.py" \
      --min_mapq 10 --max_mismatches 2 \
      "$RESULTS/tmp/$run.namesort.bam" "$RESULTS/tmp/$run.qcmarked.bam"
    samtools view -@ "$THREADS" -u -F 512 "$RESULTS/tmp/$run.qcmarked.bam" \
    | samtools sort -@ "$THREADS" -m 1G -T "$tmp_prefix.coordsort" \
        -o "$out.part" -
    mv "$out.part" "$out"
    samtools index -@ "$THREADS" "$out"
    rm -f "$RESULTS/tmp/$run.namesort.bam" "$RESULTS/tmp/$run.qcmarked.bam"
  done
fi

if run_stage call; then
  cut -f1 "$SAMPLES" | tail -n +2 | sort -u | while read -r donor; do
    mapfile -t bams < <(awk -F '\t' -v d="$donor" 'NR>1 && $1==d {print $3}' "$SAMPLES" \
      | sed "s#^#$RESULTS/bam/#; s#\$#.filtered.bam#")
    for bam in "${bams[@]}"; do
      [[ -s "$bam" && -s "$bam.bai" ]] || { echo "Missing BAM or index: $bam" >&2; exit 1; }
    done

    het="$RESULTS/variants/$donor.heterozygous.vcf.gz"
    if [[ -s "$het" && -s "$het.tbi" ]]; then
      echo "Variant calling already complete: $donor"
      continue
    fi
    echo "Calling variants for $donor from ${#bams[@]} libraries"
    merged="$RESULTS/bam/$donor.merged.bam"
    raw="$RESULTS/variants/$donor.raw.bcf"
    norm="$RESULTS/variants/$donor.norm.bcf"
    filt="$RESULTS/variants/$donor.biallelic_snps.bcf"
    samtools merge -@ "$THREADS" -f "$merged.part" "${bams[@]}"
    mv "$merged.part" "$merged"
    samtools index -@ "$THREADS" "$merged"

    bcftools mpileup --threads "$THREADS" --redo-BAQ --adjust-MQ 50 \
      --gap-frac 0.05 --max-depth 10000 -a FORMAT/AD,FORMAT/DP \
      -Ou -f "$FA" "$merged" \
    | bcftools call --threads "$THREADS" --keep-alts --multiallelic-caller \
        --format-fields GQ -Ob -o "$raw"
    bcftools norm --threads "$THREADS" --check-ref x -m - -Ob \
      -o "$norm" -f "$FA" "$raw"
    bcftools filter --threads "$THREADS" \
      -i 'QUAL>=10 && FORMAT/GQ>=20 && FORMAT/DP>=10' \
      --SnpGap 3 --IndelGap 10 -Ou "$norm" \
    | bcftools view --threads "$THREADS" -m2 -M2 -v snps -Ob -o "$filt"
    bcftools index -f "$filt"
    bcftools view -i 'GT="het" && MIN(FORMAT/AD)>=5' -Oz -o "$het.part" "$filt"
    mv "$het.part" "$het"
    bcftools index -f -t "$het"
  done
fi

if run_stage count; then
  awk -F '\t' 'NR>1 {print $1 "\t" $3}' "$SAMPLES" \
  | while IFS=$'\t' read -r donor run; do
    out="$RESULTS/counts/$run.mapped"
    if [[ -s "$out" ]]; then
      echo "Allele counts already complete: $run"
      continue
    fi
    python3 "$REPO_DIR/pipeline/count_alleles.py" \
      --bam "$RESULTS/bam/$run.filtered.bam" \
      --vcf "$RESULTS/variants/$donor.heterozygous.vcf.gz" \
      --format mixalime --sample-id "$run" \
      --output "$out.part"
    mv "$out.part" "$out"
  done
fi

echo "Stage '$STAGE' complete: $RESULTS"

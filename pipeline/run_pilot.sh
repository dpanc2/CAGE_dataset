#!/usr/bin/env bash
set -euo pipefail

# Reproduce the non-WASP primary CAGE/ASV preprocessing on the two-run pilot.
# Large files default to /mnt/newdata; override with CAGE_DATA_ROOT if needed.

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${CAGE_DATA_ROOT:-/mnt/newdata/dpanc/CAGE_dataset_data}"
THREADS="${THREADS:-16}"
SAMPLES="${SAMPLES:-$REPO_DIR/pipeline/pilot_samples.tsv}"
FASTQ_DIR="$DATA_ROOT/fastq"
REF_DIR="$DATA_ROOT/reference/gencode_v39"
RESULTS="$DATA_ROOT/results/pilot"
FA="$REF_DIR/GRCh38.primary_assembly.genome.fa"
GTF="$REF_DIR/gencode.v39.annotation.gtf"
INDEX="$REF_DIR/hisat2/GRCh38_gencode_v39"

need() { command -v "$1" >/dev/null || { echo "Missing required command: $1" >&2; exit 1; }; }
for cmd in curl gzip hisat2 hisat2-build hisat2_extract_splice_sites.py \
           hisat2_extract_exons.py samtools bcftools python3; do need "$cmd"; done

mkdir -p "$FASTQ_DIR" "$REF_DIR/hisat2" "$RESULTS"/{bam,logs,variants,counts}

if [[ ! -s "$FA" ]]; then
  curl -L --fail --retry 10 -o "$FA.gz.part" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz
  gzip -t "$FA.gz.part"
  mv "$FA.gz.part" "$FA.gz"
  gzip -dk "$FA.gz"
fi
if [[ ! -s "$GTF" ]]; then
  curl -L --fail --retry 10 -o "$GTF.gz.part" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
  gzip -t "$GTF.gz.part"
  mv "$GTF.gz.part" "$GTF.gz"
  gzip -dk "$GTF.gz"
fi

samtools faidx "$FA"
if [[ ! -s "$INDEX.1.ht2" && ! -s "$INDEX.1.ht2l" ]]; then
  hisat2_extract_splice_sites.py "$GTF" > "$REF_DIR/splice_sites.tsv"
  hisat2_extract_exons.py "$GTF" > "$REF_DIR/exons.tsv"
  hisat2-build -p "$THREADS" --ss "$REF_DIR/splice_sites.tsv" \
    --exon "$REF_DIR/exons.tsv" "$FA" "$INDEX"
fi

tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r donor run chamber fastq_name; do
  fq="$FASTQ_DIR/$fastq_name"
  out="$RESULTS/bam/$run.filtered.bam"
  [[ -s "$fq" ]] || { echo "Missing FASTQ: $fq" >&2; exit 1; }
  [[ -s "$out" ]] && continue
  hisat2 -p "$THREADS" --very-sensitive -x "$INDEX" -U "$fq" \
    --rg-id "$run" --rg "SM:$donor" \
    --summary-file "$RESULTS/logs/$run.hisat2.txt" 2> "$RESULTS/logs/$run.stderr.txt" \
  | samtools view -@ "$THREADS" -u - \
  | samtools sort -@ "$THREADS" -n -o "$RESULTS/bam/$run.namesort.bam" -

  python3 "$REPO_DIR/pipeline/vendor/filter_reads.py" \
    --min_mapq 10 --max_mismatches 2 \
    "$RESULTS/bam/$run.namesort.bam" "$RESULTS/bam/$run.qcmarked.bam"
  samtools view -@ "$THREADS" -u -F 512 "$RESULTS/bam/$run.qcmarked.bam" \
  | samtools sort -@ "$THREADS" -o "$out" -
  samtools index -@ "$THREADS" "$out"
  rm -f "$RESULTS/bam/$run.namesort.bam" "$RESULTS/bam/$run.qcmarked.bam"
done

cut -f1 "$SAMPLES" | tail -n +2 | sort -u | while read -r donor; do
  mapfile -t bams < <(awk -F '\t' -v d="$donor" 'NR>1 && $1==d {print $2}' "$SAMPLES" \
    | sed "s#^#$RESULTS/bam/#; s#\$#.filtered.bam#")
  merged="$RESULTS/bam/$donor.merged.bam"
  samtools merge -@ "$THREADS" -f "$merged" "${bams[@]}"
  samtools index -@ "$THREADS" "$merged"

  raw="$RESULTS/variants/$donor.raw.bcf"
  norm="$RESULTS/variants/$donor.norm.bcf"
  filt="$RESULTS/variants/$donor.biallelic_snps.bcf"
  het="$RESULTS/variants/$donor.heterozygous.vcf.gz"
  bcftools mpileup --threads "$THREADS" --redo-BAQ --adjust-MQ 50 \
    --gap-frac 0.05 --max-depth 10000 -a FORMAT/AD,FORMAT/DP -Ou -f "$FA" "$merged" \
  | bcftools call --threads "$THREADS" --keep-alts --multiallelic-caller -Ob -o "$raw"
  bcftools norm --threads "$THREADS" --check-ref x -m - -Ob -o "$norm" -f "$FA" "$raw"
  bcftools filter --threads "$THREADS" -i 'QUAL>=10 && FORMAT/GQ>=20 && FORMAT/DP>=10' \
    --SnpGap 3 --IndelGap 10 -Ou "$norm" \
  | bcftools view --threads "$THREADS" -m2 -M2 -v snps -Ob -o "$filt"
  bcftools index -f "$filt"
  bcftools view -i 'GT="het" && MIN(FORMAT/AD)>=5' -Oz -o "$het" "$filt"
  bcftools index -f -t "$het"

  awk -F '\t' -v d="$donor" 'NR>1 && $1==d {print $2}' "$SAMPLES" | while read -r run; do
    python3 "$REPO_DIR/pipeline/count_alleles.py" \
      --bam "$RESULTS/bam/$run.filtered.bam" --vcf "$het" \
      --format mixalime --sample-id "$run" \
      --output "$RESULTS/counts/$run.mapped"
  done
done

echo "Pilot complete: $RESULTS"

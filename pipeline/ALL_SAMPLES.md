# All 109 CAGE libraries

This workflow downloads and preprocesses all GSE150736 libraries and uses the
published subject numbers from Deviatiiarov et al. Supplementary Table 1 to
merge libraries from the same heart.

On `airigpu`:

```bash
cd /home/dpanc/CAGE_dataset
conda activate cage-asv

mkdir -p scripts/logs
nohup env THREADS=8 STAGE=all bash pipeline/run_all.sh \
  > scripts/logs/run_all.log 2>&1 &
echo $! > scripts/logs/run_all.pid
```

Follow progress with:

```bash
tail -f scripts/logs/run_all.log
```

The default data root is `/mnt/newdata/dpanc/CAGE_dataset_data`. Override it
with `CAGE_DATA_ROOT` if necessary.

The script is restartable. Completed FASTQs, filtered BAMs, donor VCFs, and
count tables are skipped. To run or resume only one stage:

```bash
THREADS=8 STAGE=download bash pipeline/run_all.sh
THREADS=8 STAGE=align    bash pipeline/run_all.sh
THREADS=8 STAGE=call     bash pipeline/run_all.sh
THREADS=8 STAGE=count    bash pipeline/run_all.sh
```

Main outputs:

```text
reference/heart_atlas_metadata/GSE150736_all_samples.tsv
results/all/bam/<run>.filtered.bam
results/all/bam/subject_<NN>.merged.bam
results/all/variants/subject_<NN>.heterozygous.vcf.gz
results/all/counts/<run>.allele_counts.tsv
```

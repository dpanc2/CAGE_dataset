# Pilot pipeline

This reproduces the paper's primary non-WASP workflow for the inferred two-run
pilot donor: GENCODE v39/hg38 HISAT2 alignment, StamLab read filtering,
donor-level bcftools calling/filtering, and sample-level REF/ALT counts.

Create the software environment and move the two FASTQs to the large data area:

```bash
./pipeline/setup_environment.sh
mkdir -p /mnt/newdata/dpanc/CAGE_dataset_data/fastq
mv data/fastq/SRR1650629{7,8}.fastq.gz /mnt/newdata/dpanc/CAGE_dataset_data/fastq/
MAMBA_ROOT_PREFIX=/mnt/newdata/dpanc/micromamba \
  /mnt/newdata/dpanc/micromamba/bin/micromamba run -n cage-asv \
  env THREADS=16 bash pipeline/run_pilot.sh
```

The dbSNP v151 annotation and optional WASP remapping are deliberately separate
follow-up stages: neither is required to obtain the initial heterozygous loci and
allelic counts, and both require additional versioned reference inputs.

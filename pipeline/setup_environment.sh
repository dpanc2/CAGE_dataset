#!/usr/bin/env bash
set -euo pipefail

ROOT="${MAMBA_ROOT_PREFIX:-/mnt/newdata/dpanc/micromamba}"
BIN="$ROOT/bin/micromamba"
mkdir -p "$ROOT/bin"
if [[ ! -x "$BIN" ]]; then
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C "$ROOT/bin" --strip-components=1 bin/micromamba
fi
export MAMBA_ROOT_PREFIX="$ROOT"
"$BIN" create -y -n cage-asv -f "$(dirname "$0")/environment.yml"
"$BIN" run -n cage-asv bash -c \
  'hisat2 --version | head -n1; samtools --version | head -n1; bcftools --version | head -n1; python -c "import pysam; print(\"pysam\", pysam.__version__)"'


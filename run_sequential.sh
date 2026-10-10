#!/usr/bin/env sh
set -eu
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
cd "$(dirname "$0")"
exec Rscript run_analysis.R "${1:-reference}"

#!/bin/bash
#SBATCH --job-name=partition
#SBATCH --output=%x_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
# Size the array with `liposome-partition count $ROOTS` (--array=0-N-1).

set -euo pipefail

if [ -n "${GMX_SETUP:-}" ]; then
    source "$GMX_SETUP"
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
ROOTS="${ROOTS:-.}"
CACHE_DIR="${CACHE_DIR:-results_cache}"

liposome-partition compute $ROOTS --index "$SLURM_ARRAY_TASK_ID" --cache-dir "$CACHE_DIR"

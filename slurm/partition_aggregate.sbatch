#!/bin/bash
#SBATCH --job-name=partition_aggregate
#SBATCH --output=%x_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb

set -euo pipefail

CACHE_DIR="${CACHE_DIR:-results_cache}"
OUTPUT_DIR="${OUTPUT_DIR:-partition_results}"

liposome-partition aggregate --cache-dir "$CACHE_DIR" --output-dir "$OUTPUT_DIR"

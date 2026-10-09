#!/bin/bash
# Parallel partition-ratio analysis: one array task per system, then one aggregation job.
#   ROOTS="run1 run2" CACHE_DIR=cache THROTTLE=10 slurm/submit_partition.sh [sbatch options...]
set -euo pipefail

ROOTS="${ROOTS:-.}"
CACHE_DIR="${CACHE_DIR:-results_cache}"
THROTTLE="${THROTTLE:-20}"
here="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$CACHE_DIR"

n=$(liposome-partition count $ROOTS)
if [ "$n" -eq 0 ]; then
    echo "No drug-loaded systems found under: $ROOTS" >&2
    exit 1
fi
echo "Found $n systems; submitting array 0-$((n - 1))%$THROTTLE"

array_id=$(ROOTS="$ROOTS" CACHE_DIR="$CACHE_DIR" \
    sbatch --parsable --array="0-$((n - 1))%$THROTTLE" "$@" "$here/partition_array.sbatch")
echo "Array job: $array_id"

agg_id=$(CACHE_DIR="$CACHE_DIR" \
    sbatch --parsable --dependency="afterok:$array_id" "$@" "$here/partition_aggregate.sbatch")
echo "Aggregation job: $agg_id"

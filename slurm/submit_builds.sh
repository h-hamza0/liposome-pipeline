#!/bin/bash
# Submit one build job per configuration file in a directory.
#   slurm/submit_builds.sh CONFIG_DIR [sbatch options...]
set -euo pipefail

dir="${1:?usage: submit_builds.sh CONFIG_DIR [sbatch options]}"
shift
here="$(cd "$(dirname "$0")" && pwd)"

for config in "$dir"/*.txt; do
    echo "submitting $config"
    sbatch "$@" "$here/build.sbatch" "$config"
done

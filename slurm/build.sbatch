#!/bin/bash
#SBATCH --job-name=liposome-build
#SBATCH --output=%x_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32gb
#SBATCH --gres=gpu:1
# Site-specific options belong on the command line:
#   GMX_SETUP=/path/to/GMXRC sbatch --partition=gpu slurm/build.sbatch config.txt

set -euo pipefail

if [ -n "${GMX_SETUP:-}" ]; then
    source "$GMX_SETUP"
fi

export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

liposome-build run "$1" --ntomp "${SLURM_CPUS_PER_TASK:-12}" --ntmpi 1 --nb gpu --pin off

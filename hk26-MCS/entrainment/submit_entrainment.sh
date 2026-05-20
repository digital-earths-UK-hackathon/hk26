#!/bin/bash
# Submit CAPE/CIN/w_norm computation for the WAM region as a SLURM array.
#
# Run --init first, then submit this script:
#   python calc_entrainment_wam.py --init
#   sbatch submit_cape_cin.sh
#
# Array range 0-40 covers 3249 timesteps at 80 per chunk (ceil(3249/80) = 41 jobs).

#SBATCH --job-name=entrainment_wam
#SBATCH --array=0-40
#SBATCH --account=hrcm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --time=10:0:00
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --output=logs/entrainment_%A_%a.out
#SBATCH --error=logs/entrainment_%A_%a.err

set -euo pipefail

mkdir -p logs

# Inherits current env.
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate hk26_env

python calc_entrainment_wam.py --chunk "$SLURM_ARRAY_TASK_ID"

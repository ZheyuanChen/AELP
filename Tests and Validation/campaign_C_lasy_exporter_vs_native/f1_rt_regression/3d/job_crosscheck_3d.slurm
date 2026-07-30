#!/usr/bin/env bash
# 3D file-level crosscheck (see crosscheck_original_files.py). Memmap
# slab scan over 4 x 5.8 GB files -- I/O-bound, small resident memory.
#SBATCH --job-name=campC3-3d-crosscheck
#SBATCH --partition=nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --account=pet-pic-2022
#SBATCH --mail-type=NONE
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

set -uo pipefail
source /mnt/scratch/users/pnd531/project_epoch/.venv/bin/activate
cd "${SLURM_SUBMIT_DIR:?Submit this job from f1_rt_regression/3d}"
echo "Host: $(hostname)   Start: $(date)"
timeout 2100 python crosscheck_original_files.py
echo "crosscheck_rc=$?"
echo "Done: $(date)"

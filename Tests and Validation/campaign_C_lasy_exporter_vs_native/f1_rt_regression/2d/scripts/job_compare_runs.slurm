#!/usr/bin/env bash
# Run-level 2D comparison (see compare_runs.py). Loads 2 x 16 small 2D SDFs.
#SBATCH --job-name=campC3-2d-compare
#SBATCH --partition=nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --account=pet-pic-2022
#SBATCH --mail-type=NONE
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

set -uo pipefail
source /mnt/scratch/users/pnd531/project_epoch/.venv/bin/activate
cd "${SLURM_SUBMIT_DIR:?Submit this job from f1_rt_regression/2d/scripts}"
echo "Host: $(hostname)   Start: $(date)"
timeout 1500 python compare_runs.py
echo "compare_rc=$?"
echo "Done: $(date)"

#!/usr/bin/env bash
# Run-level 3D comparison (see compare_runs_3d.py). One 5.9 GB SDF at a
# time from each of 2 x 16 snapshots -- I/O heavy, moderate memory.
#SBATCH --job-name=campC3-3d-compare
#SBATCH --partition=nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --account=pet-pic-2022
#SBATCH --mail-type=NONE
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

set -uo pipefail
source /mnt/scratch/users/pnd531/project_epoch/.venv/bin/activate
cd "${SLURM_SUBMIT_DIR:?Submit this job from f1_rt_regression/3d/scripts}"
echo "Host: $(hostname)   Start: $(date)"
timeout 3000 python compare_runs_3d.py
echo "compare_rc=$?"
echo "Done: $(date)"

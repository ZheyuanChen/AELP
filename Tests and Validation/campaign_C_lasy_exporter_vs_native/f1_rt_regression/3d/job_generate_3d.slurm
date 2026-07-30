#!/usr/bin/env bash
# 3D rt-exporter generation (see generate_exporter_rt.py). ~25 GB peak
# RAM (two 5.8 GB output arrays + transpose copies), minutes of numpy.
#SBATCH --job-name=campC3-3d-generate
#SBATCH --partition=nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:30:00
#SBATCH --account=pet-pic-2022
#SBATCH --mail-type=NONE
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

set -uo pipefail
source /mnt/scratch/users/pnd531/project_epoch/.venv/bin/activate
cd "${SLURM_SUBMIT_DIR:?Submit this job from f1_rt_regression/3d}"
echo "Host: $(hostname)   Start: $(date)"
timeout 4800 python generate_exporter_rt.py
echo "generate_rc=$?"
echo "Done: $(date)"

#!/bin/bash
#SBATCH --job-name=mutect
#SBATCH --output=mutect_%j.out
#SBATCH --error=mutect_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH -p gpu
# SBATCH -A YOUR_ACCOUNT   # <— leave commented in public template

module load apptainer
module load nextflow
export NXF_SINGULARITY_CMD=apptainer
# Optional binds via env (prefer this over hard-coding in config)
export SINGULARITY_RUN_OPTS="--nv"

nextflow run . -profile slurm -resume \
  --samplesheet path/to/sample.csv \
  --outdir results/run1 \
  --num_gpus 1 \
  --rg_lb sureselect_v8 --rg_pl ILLUMINA
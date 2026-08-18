#!/bin/bash
#SBATCH --job-name=somatic_nf
#SBATCH --output=somatic_nf_%j.out
#SBATCH --error=somatic_nf_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -p compute
#SBATCH -A o240001

source ~/.bashrc
ml apptainer
ml nextflow
export NXF_SINGULARITY_CMD=apptainer

# Defaults now use the Broad hg38/assembly38 reference/resource bundle in main.nf.
# Override via REF_FA / DBSNP_VCF / MILLS_VCF / MUTECT_GERMLINE_RESOURCE_VCF /
# CONTAMINATION_SITES_VCF / MUTECT_PON_VCF if you need a different bundle.
nextflow run main.nf \
  -profile slurm \
  -resume \
  --samplesheet sample_cram_copy.csv \
  --slurm_account o250039 \
  --outdir /project/o240001_SBUFF67/ALL_FF68/exome/rerun_cram_2 \
  --num_gpus 1
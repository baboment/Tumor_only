#!/bin/bash
#SBATCH --job-name=somatic_nf_hg19
#SBATCH --output=somatic_nf_hg19_%j.out
#SBATCH --error=somatic_nf_hg19_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -p compute
#SBATCH -A o250038

source ~/.bashrc
ml apptainer
ml nextflow
export NXF_SINGULARITY_CMD=apptainer

# hg19/b37 run: SP39 (proband), SP40 (father), SP41 (mother) - trio exome,
# three independent tumour-only rows, no pairing.
#
# -profile slurm,hg19 supplies every reference path from conf/hg19.config.
# Do NOT also export REF_FA / DBSNP_VCF / MUTECT_PON_VCF here: the profile
# overrides them, so they would be silently ignored.
#
# Prerequisites, once per cluster:
#   sbatch scripts/fetch_somatic_b37.sh     # b37 Mutect2 + ClairS-TO resources
#   ./scripts/prepare_bed_b37.sh            # S31285117 hg19 BEDs -> b37
#
# --slurm_account is the billed account. A bare --account is silently ignored
# by this pipeline (there is no params.account); keep it and #SBATCH -A in sync.

nextflow run main.nf \
  -profile slurm,hg19 \
  -resume \
  --slurm_account o250038 \
  --samplesheet /project/o240001_SBUFF67/script_zam/Tumor_only/sample_trio_b37.csv \
  --outdir /project/o250038_ALL_MET/case-report/tumor_only_hg19 \
  --num_gpus 1 \
  --rg_lb sureselect_v7 \
  --rg_pl ILLUMINA

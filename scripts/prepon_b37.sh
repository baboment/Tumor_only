#!/bin/bash
#SBATCH --job-name=prepon_b37
#SBATCH --output=prepon_b37_%j.out
#SBATCH --error=prepon_b37_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH -p gpu
#SBATCH -A o250038
# =============================================================================
# Build the Parabricks PON index for the b37 Mutect2 panel of normals.
#
#   sbatch scripts/prepon_b37.sh
#
# WHY THIS IS A SEPARATE STEP. `pbrun mutectcaller --pon x.vcf.gz` does not read
# the VCF directly - it requires a sidecar `x.vcf.gz.pon` produced by
# `pbrun prepon`, and dies with
#     <pon>.pon not found. Make sure you run prepon first
# if it is absent. The hg38 PONs in /common/db/human_ref/hg38/somatic_hg38 each
# carry one; the b37 panel fetched by scripts/fetch_somatic_b37.sh does not,
# because prepon needs a GPU and that script runs on the compute partition.
#
# prepon runs in place, writing <pon>.pon beside the input, so conf/hg19.config
# needs no change. Unlike scripts/prepare_mutect_pon.sh there is no reheader or
# contig-subset step here: fetch_somatic_b37.sh already reheadered this panel
# against the reference .fai and asserted its contigs.
# =============================================================================
set -euo pipefail

PON_VCF=${PON_VCF:-/common/db/human_ref/hg19/somatic_b37/Mutect2-exome-panel.b37.vcf.gz}
PB_SIF=${PB_SIF:-/common/sif/clara-parabricks/4.5.1.sif}
NUM_GPUS=${NUM_GPUS:-1}
APPTAINER_BIND_PATHS=${APPTAINER_BIND_PATHS:-/common:/common,/project:/project}

[[ -f "$PON_VCF" ]] || { echo "PON not found: $PON_VCF" >&2; exit 1; }
[[ -f "$PB_SIF"  ]] || { echo "Parabricks container not found: $PB_SIF" >&2; exit 1; }

if [[ -s "$PON_VCF.pon" ]]; then
  echo "[skip ] $PON_VCF.pon already exists"
  exit 0
fi

module load apptainer 2>/dev/null || true
command -v apptainer >/dev/null || { echo "apptainer not on PATH" >&2; exit 1; }

binds=()
for p in ${APPTAINER_BIND_PATHS//,/ }; do binds+=(--bind "$p"); done

echo "[prepon] $PON_VCF"
apptainer exec --nv --cleanenv "${binds[@]}" "$PB_SIF" \
  pbrun prepon --num-gpus "$NUM_GPUS" --in-pon-file "$PON_VCF"

[[ -s "$PON_VCF.pon" ]] || { echo "prepon finished but $PON_VCF.pon was not created" >&2; exit 1; }
ls -la "$PON_VCF" "$PON_VCF.pon"
echo "PON index ready."

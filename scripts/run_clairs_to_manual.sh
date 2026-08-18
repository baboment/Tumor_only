#!/bin/bash
#SBATCH --job-name=clairs_to_manual
#SBATCH --output=clairs_to_manual_%j.out
#SBATCH --error=clairs_to_manual_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=512G
#SBATCH --cpus-per-task=48
#SBATCH -p compute
#SBATCH -A o240001

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_clairs_to_manual.sh [options] [-- extra ClairS-TO args]
  sbatch scripts/run_clairs_to_manual.sh [options] [-- extra ClairS-TO args]

Default values are filled in for the current Tumor_only run.

Options:
  --sample NAME          Sample name. Default: DC301
  --cram PATH            Input tumor CRAM. Default: <outroot>/cram_bqsr/<sample>.bqsr.cram
  --crai PATH            Input CRAI. Default: <outroot>/cram_bqsr/<sample>.bqsr.cram.crai
  --outdir PATH          Output directory. Default: <outroot>/manual_clairs_to/<sample>
  --ref PATH             Reference FASTA.
  --bed PATH             BED file for ROI calling.
  --platform NAME        ClairS-TO platform. Default: ilmn_ssrs
  --threads N            Threads for run_clairs_to. Default: 48
  --container PATH       ClairS-TO container path.
  --fast                 Add --disable_intermediate_phasing.
  --debug                Add --disable_intermediate_phasing and --disable_verdict.
  --dry-run              Print commands and exit.
  --help                 Show this help.

Anything after `--` is passed directly to run_clairs_to.
EOF
}

timestamp() {
  date '+%F %T'
}

SAMPLE="DC301"
OUTROOT="/project/o240001_SBUFF67/ALL_FF67/normal/tumor_only"
REF="/common/db/human_ref/hg38/v0/Homo_sapiens_assembly38.fasta"
BED="/scratch/o240001_SBUFF67/exome/03_preprocess_result/S33266436_Padded.chr.clean.bed"
PLATFORM="ilmn_ssrs"
THREADS="48"
CONTAINER="/common/sif/Clairs-TO/Clairs-TO.v0.4.2.sif"
PON_GNOMAD="/project/o240001_SBUFF67/ClairS-TO/PON/gnomad.r2.1.af-ge-0.001.sites.vcf.gz"
PON_DBSNP="/project/o240001_SBUFF67/ClairS-TO/PON/dbsnp.b138.non-somatic.sites.vcf.gz"
PON_1000G="/project/o240001_SBUFF67/ClairS-TO/PON/1000g-pon.sites.vcf.gz"
PON_COLORSDB="/project/o240001_SBUFF67/ClairS-TO/PON/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz"
FAST_MODE=0
DEBUG_MODE=0
DRY_RUN=0
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample)
      SAMPLE="$2"
      shift 2
      ;;
    --cram)
      CRAM="$2"
      shift 2
      ;;
    --crai)
      CRAI="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --ref)
      REF="$2"
      shift 2
      ;;
    --bed)
      BED="$2"
      shift 2
      ;;
    --platform)
      PLATFORM="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --container)
      CONTAINER="$2"
      shift 2
      ;;
    --fast)
      FAST_MODE=1
      shift
      ;;
    --debug)
      DEBUG_MODE=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    --)
      shift
      EXTRA_ARGS+=("$@")
      break
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

CRAM="${CRAM:-${OUTROOT}/cram_bqsr/${SAMPLE}.bqsr.cram}"
CRAI="${CRAI:-${OUTROOT}/cram_bqsr/${SAMPLE}.bqsr.cram.crai}"
OUTDIR="${OUTDIR:-${OUTROOT}/manual_clairs_to/${SAMPLE}}"
LOGDIR="${OUTDIR}/logs"

[[ -f "$CRAM" ]] || { echo "Missing CRAM: $CRAM" >&2; exit 1; }
[[ -f "$CRAI" ]] || { echo "Missing CRAI: $CRAI" >&2; exit 1; }
[[ -f "$REF" ]] || { echo "Missing reference FASTA: $REF" >&2; exit 1; }
[[ -f "$BED" ]] || { echo "Missing BED file: $BED" >&2; exit 1; }
[[ -f "$CONTAINER" ]] || { echo "Missing container: $CONTAINER" >&2; exit 1; }
[[ -f "$PON_GNOMAD" ]] || { echo "Missing PoN: $PON_GNOMAD" >&2; exit 1; }
[[ -f "$PON_DBSNP" ]] || { echo "Missing PoN: $PON_DBSNP" >&2; exit 1; }
[[ -f "$PON_1000G" ]] || { echo "Missing PoN: $PON_1000G" >&2; exit 1; }
[[ -f "$PON_COLORSDB" ]] || { echo "Missing PoN: $PON_COLORSDB" >&2; exit 1; }

set +u
source ~/.bashrc
set -u
module load apptainer

mkdir -p "$OUTDIR" "$LOGDIR"

BIND_ARGS=(
  --bind /common:/common
  --bind /project:/project
  --bind /scratch:/scratch
)

RUN_ARGS=(
  --tumor_bam_fn "$CRAM"
  --ref_fn "$REF"
  --sample_name "$SAMPLE"
  --threads "$THREADS"
  --platform "$PLATFORM"
  --panel_of_normals "${PON_GNOMAD},${PON_DBSNP},${PON_1000G},${PON_COLORSDB}"
  --panel_of_normals_require_allele_matching "True,True,False,False"
  --output_dir "$OUTDIR"
  --conda_prefix /opt/micromamba/envs/clairs-to
  --bed_fn "$BED"
)

if (( FAST_MODE )); then
  RUN_ARGS+=(--disable_intermediate_phasing)
fi

if (( DEBUG_MODE )); then
  RUN_ARGS+=(--disable_intermediate_phasing --disable_verdict)
fi

if (( ${#EXTRA_ARGS[@]} > 0 )); then
  RUN_ARGS+=("${EXTRA_ARGS[@]}")
fi

RUN_CMD=(
  apptainer exec
  --cleanenv
  "${BIND_ARGS[@]}"
  "$CONTAINER"
  run_clairs_to
  "${RUN_ARGS[@]}"
)

MERGE_CMD=(
  apptainer exec
  --cleanenv
  "${BIND_ARGS[@]}"
  "$CONTAINER"
  bcftools concat -a "${OUTDIR}/snv.vcf.gz" "${OUTDIR}/indel.vcf.gz"
)

SORT_CMD=(
  apptainer exec
  --cleanenv
  "${BIND_ARGS[@]}"
  "$CONTAINER"
  bcftools sort -Oz -o "${OUTDIR}/${SAMPLE}.clairs.vcf.gz"
)

INDEX_CMD=(
  apptainer exec
  --cleanenv
  "${BIND_ARGS[@]}"
  "$CONTAINER"
  bcftools index -t "${OUTDIR}/${SAMPLE}.clairs.vcf.gz"
)

{
  echo "# Generated at $(timestamp)"
  echo "# Sample: ${SAMPLE}"
  echo "# Output: ${OUTDIR}"
  printf '%q ' "${RUN_CMD[@]}"
  echo
} > "${LOGDIR}/run_clairs_to.command.sh"

echo "[$(timestamp)] Sample: ${SAMPLE}"
echo "[$(timestamp)] Output dir: ${OUTDIR}"

if (( DRY_RUN )); then
  cat "${LOGDIR}/run_clairs_to.command.sh"
  exit 0
fi

SECONDS=0
echo "[$(timestamp)] Starting ClairS-TO"
"${RUN_CMD[@]}" 2>&1 | tee "${LOGDIR}/run_clairs_to.log"

SNV_VCF="${OUTDIR}/snv.vcf.gz"
INDEL_VCF="${OUTDIR}/indel.vcf.gz"
FINAL_VCF="${OUTDIR}/${SAMPLE}.clairs.vcf.gz"

if [[ -f "$SNV_VCF" && -f "$INDEL_VCF" ]]; then
  echo "[$(timestamp)] Merging SNV and indel VCFs"
  "${MERGE_CMD[@]}" | "${SORT_CMD[@]}"
elif [[ -f "$SNV_VCF" ]]; then
  echo "[$(timestamp)] Only SNV VCF present; copying to final output"
  cp -f "$SNV_VCF" "$FINAL_VCF"
elif [[ -f "$INDEL_VCF" ]]; then
  echo "[$(timestamp)] Only indel VCF present; copying to final output"
  cp -f "$INDEL_VCF" "$FINAL_VCF"
else
  echo "ClairS-TO finished without producing snv.vcf.gz or indel.vcf.gz in ${OUTDIR}" >&2
  exit 1
fi

echo "[$(timestamp)] Indexing final VCF"
"${INDEX_CMD[@]}"

echo "[$(timestamp)] Finished in ${SECONDS}s"
echo "[$(timestamp)] Final VCF: ${FINAL_VCF}"

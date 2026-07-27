#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Prepare a Mutect PON against the exact reference .fai used by this repo.

This script does not modify the pipeline. It:
1. makes a new reheadered copy of the input PON with bcftools reheader --fai
2. tabix-indexes the reheadered copy
3. reruns pbrun prepon on that copy

Usage:
  bash scripts/prepare_mutect_pon.sh \
    --in-pon-file /path/to/input.pon.vcf.gz \
    --ref-fai /path/to/reference.fa.fai \
    --out-pon-file /path/to/input.pon.reheadered.vcf.gz

Optional:
  --pb-container /path/to/clara-parabricks.sif
  --num-gpus 1

Environment fallbacks:
  PB_SIF    Parabricks container path when pbrun is not on PATH
  NUM_GPUS  GPU count for prepon (default: 1)
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

cleanup() {
  if [[ -n "${tmpdir:-}" && -d "${tmpdir}" ]]; then
    rm -rf "${tmpdir}"
  fi
}

maybe_source_bashrc() {
  if [[ -f "${HOME}/.bashrc" ]]; then
    local had_nounset=0
    if [[ $- == *u* ]]; then
      had_nounset=1
      set +u
    fi

    # shellcheck disable=SC1090
    source "${HOME}/.bashrc"

    if [[ ${had_nounset} -eq 1 ]]; then
      set -u
    fi
  fi
}

maybe_load_module() {
  local cmd="$1"
  local module_name="$2"

  if command -v "${cmd}" >/dev/null 2>&1; then
    return 0
  fi

  if ! command -v module >/dev/null 2>&1 && ! command -v ml >/dev/null 2>&1; then
    maybe_source_bashrc
  fi

  if command -v module >/dev/null 2>&1; then
    module load "${module_name}" >/dev/null 2>&1 || true
  elif command -v ml >/dev/null 2>&1; then
    ml "${module_name}" >/dev/null 2>&1 || true
  fi
}

canonical_dir() {
  local path="$1"
  local dir

  dir="$(dirname "${path}")"
  (
    cd "${dir}" >/dev/null 2>&1
    pwd -P
  )
}

canonical_path() {
  local path="$1"
  local dir
  local base

  dir="$(canonical_dir "${path}")"
  base="$(basename "${path}")"
  echo "${dir}/${base}"
}

run_prepon() {
  local pon_vcf="$1"
  local num_gpus="$2"
  local pb_container="${3:-}"

  if command -v pbrun >/dev/null 2>&1; then
    pbrun prepon \
      --num-gpus "${num_gpus}" \
      --in-pon-file "${pon_vcf}"
    return 0
  fi

  [[ -n "${pb_container}" ]] || die "pbrun is not on PATH and no --pb-container/PB_SIF was provided."

  maybe_load_module apptainer apptainer
  command -v apptainer >/dev/null 2>&1 || die "apptainer is required when pbrun is not on PATH."
  [[ -f "${pb_container}" ]] || die "Parabricks container not found: ${pb_container}"

  local bind_dir
  bind_dir="$(canonical_dir "${pon_vcf}")"

  apptainer exec \
    --nv \
    --bind "${bind_dir}:${bind_dir}" \
    "${pb_container}" \
    pbrun prepon \
      --num-gpus "${num_gpus}" \
      --in-pon-file "${pon_vcf}"
}

input_pon=''
ref_fai=''
output_pon=''
pb_container="${PB_SIF:-}"
num_gpus="${NUM_GPUS:-1}"
tmpdir=''

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in-pon-file)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      input_pon="$2"
      shift 2
      ;;
    --ref-fai)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      ref_fai="$2"
      shift 2
      ;;
    --out-pon-file)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      output_pon="$2"
      shift 2
      ;;
    --pb-container)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      pb_container="$2"
      shift 2
      ;;
    --num-gpus)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      num_gpus="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

[[ -n "${input_pon}" ]] || die "Missing --in-pon-file"
[[ -n "${ref_fai}" ]] || die "Missing --ref-fai"
[[ -f "${input_pon}" ]] || die "Input PON not found: ${input_pon}"
[[ -f "${ref_fai}" ]] || die "Reference .fai not found: ${ref_fai}"

input_pon="$(canonical_path "${input_pon}")"
ref_fai="$(canonical_path "${ref_fai}")"

if [[ -z "${output_pon}" ]]; then
  input_base="$(basename "${input_pon}")"
  [[ "${input_base}" == *.vcf.gz ]] || die "Input PON must end with .vcf.gz"
  output_pon="$(dirname "${input_pon}")/${input_base%.vcf.gz}.reheadered.vcf.gz"
fi

[[ "${output_pon}" != "${input_pon}" ]] || die "Output PON must be a new file, not the input file."

maybe_load_module bcftools bcftools
command -v bcftools >/dev/null 2>&1 || die "bcftools is required."

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/prepare_mutect_pon.XXXXXX")"
trap cleanup EXIT

ref_regions="${tmpdir}/ref.regions.tsv"
subset_pon="${tmpdir}/subset.input.vcf.gz"

awk 'BEGIN{OFS="\t"} {print $1,1,$2}' "${ref_fai}" > "${ref_regions}"

echo "Subsetting ${input_pon} to contigs present in ${ref_fai}"
bcftools view \
  --regions-file "${ref_regions}" \
  --output-type z \
  --output-file "${subset_pon}" \
  "${input_pon}"

bcftools index --tbi --force "${subset_pon}"

subset_records="$(bcftools index --nrecords "${subset_pon}")"
if [[ "${subset_records}" == "0" ]]; then
  die "No PON records remained after subsetting to reference contigs. The PON likely needs chromosome renaming, not just reheadering."
fi

mkdir -p "$(dirname "${output_pon}")"
output_pon="$(canonical_path "${output_pon}")"

echo "Reheadering ${subset_pon}"
bcftools reheader --fai "${ref_fai}" "${subset_pon}" > "${output_pon}"

echo "Indexing ${output_pon}"
bcftools index --tbi --force "${output_pon}"

echo "Running prepon on ${output_pon}"
run_prepon "${output_pon}" "${num_gpus}" "${pb_container}"

cat <<EOF

Prepared Mutect PON:
  VCF: ${output_pon}
  TBI: ${output_pon}.tbi
  PON index: ${output_pon}.pon

Use it by passing:
  --mutect_pon_vcf ${output_pon}
EOF

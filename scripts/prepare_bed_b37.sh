#!/bin/bash
# =============================================================================
# Convert the Agilent SureSelect V7 (S31285117) capture BEDs from UCSC hg19
# naming to b37, for the `hg19` profile.
#
#   ./scripts/prepare_bed_b37.sh
#
# The vendor bundle in /project/o250038_ALL_MET/S31285117_hg19 is UCSC-style:
# two `browser`/`track` header lines per file and chr-prefixed contigs
# (chr1..chr22, chrX, chrY - no random/hap/M). Our BAMs and every b37 resource
# use bare 1..22,X,Y,MT. hg19 and b37 share coordinates on these contigs, so the
# conversion is a rename, not a liftover.
#
# This is the b37 counterpart of
#     /project/o240001_SBUFF67/ALL_FF68/exome/code/04_prepare_bed.sh
# with that script's two hg38-isms inverted: it force-ADDS a chr prefix and
# derives .keep_chroms with `grep -E '^chr([0-9]+|X|Y|M)$'`. Here we strip the
# prefix and keep 1..22,X,Y,MT. Its GENCODE CDS section is dropped - that path
# is hardcoded to an hg38 GTF and is not needed for the capture intervals.
#
# Per the Agilent README and the convention nextflow.config enforces:
#     Regions.bed -> TARGET (params.exome_bed)      depth denominator
#     Covered.bed -> BAIT   (params.exome_bait_bed) probes + >=95% homology
# bait == target silently corrupts CollectHsMetrics enrichment numbers.
# =============================================================================
set -euo pipefail

DESIGN_DIR=${DESIGN_DIR:-/project/o250038_ALL_MET/S31285117_hg19}
BED_DIR=${BED_DIR:-$DESIGN_DIR/b37}
REF=${REF:-/common/db/human_ref/hg19/v0/Homo_sapiens_assembly19.fasta}
REF_DICT=${REF_DICT:-/common/db/human_ref/hg19/v0/Homo_sapiens_assembly19.dict}
REF_LABEL=${REF_LABEL:-b37_Broad_assembly19}
TARGET_BED=${TARGET_BED:-$DESIGN_DIR/S31285117_Regions.bed}
BAIT_BED=${BAIT_BED:-$DESIGN_DIR/S31285117_Covered.bed}
FAI="$REF.fai"

GATK_SIF=${GATK_SIF:-/common/sif/gatk/4.6.1.0.sif}
APPTAINER_BIND_PATHS=${APPTAINER_BIND_PATHS:-/common:/common,/project:/project}

module load bedtools/2.31.1 apptainer 2>/dev/null || true
command -v bedtools  >/dev/null || { echo "bedtools not on PATH"  >&2; exit 1; }
command -v apptainer >/dev/null || { echo "apptainer not on PATH" >&2; exit 1; }
[[ -f "$GATK_SIF" ]] || { echo "GATK container not found: $GATK_SIF" >&2; exit 1; }

# GATK via the container, matching TARGET_COVERAGE_METRICS in main.nf. The
# `gatk` module's launcher shells out to `python`, which does not exist on this
# cluster (only python3), so the module route fails before GATK even starts.
gatk_run() {
  local binds=()
  local p
  for p in ${APPTAINER_BIND_PATHS//,/ }; do binds+=(--bind "$p"); done
  apptainer exec --cleanenv "${binds[@]}" "$GATK_SIF" gatk "$@"
}

mkdir -p "$BED_DIR"

# b37 primary contigs. The hg38 script's regex was '^chr([0-9]+|X|Y|M)$'.
cut -f1 "$FAI" | grep -E '^([0-9]+|X|Y|MT)$' | sort -u > "$BED_DIR/.keep_chroms"

# ---------------------------------------------------------------------------
# norm_bed <in> <out>
#   strip track/browser lines -> strip chr prefix (chrM -> MT) -> keep primary
#   contigs -> sort -> merge overlapping intervals
#
# The contig filter is an exact match on column 1 inside awk, NOT the hg38
# script's `grep -wFf .keep_chroms`. With bare b37 names that grep would match
# a chromosome name appearing anywhere on the line - "1" as a word occurs in
# plenty of coordinate fields - and would let foreign contigs through.
# ---------------------------------------------------------------------------
norm_bed() {
  awk -v OFS='\t' -v keep="$BED_DIR/.keep_chroms" '
    BEGIN { while ((getline c < keep) > 0) if (c != "") ok[c] = 1 }
    /^(track|browser|#|@)/ { next }
    NF >= 3 && $2 ~ /^[0-9]+$/ {
      c = $1
      sub(/^chr/, "", c)
      if (c == "M") c = "MT"
      if (c in ok) print c, $2, $3
    }' "$1" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > "$2"
}

norm_bed "$TARGET_BED" "$BED_DIR/target.norm.b37.bed"
norm_bed "$BAIT_BED"   "$BED_DIR/bait.norm.b37.bed"

for n in target bait; do
  [[ -s "$BED_DIR/${n}.norm.b37.bed" ]] || { echo "ERROR: ${n}.norm.b37.bed is empty" >&2; exit 1; }
done

if cmp -s "$BED_DIR/target.norm.b37.bed" "$BED_DIR/bait.norm.b37.bed"; then
  echo "ERROR: target and bait BEDs are identical - CollectHsMetrics enrichment would be meaningless" >&2
  exit 1
fi

# BedToIntervalList validates every contig against the sequence dictionary and
# hard-fails on an unknown one, so this doubles as the conversion's proof.
for n in target bait; do
  gatk_run --java-options -Xmx8g BedToIntervalList \
    -I "$BED_DIR/${n}.norm.b37.bed" \
    -SD "$REF_DICT" \
    -O "$BED_DIR/${n}.b37.interval_list" \
    --QUIET true >/dev/null
done

{
  echo "Capture intervals prepared $(date -Is)"
  echo "Design    : Agilent SureSelect Human All Exon V7 (S31285117), converted hg19 -> b37"
  echo "Source    : $DESIGN_DIR"
  echo "Reference : $REF_LABEL  ($REF)"
  printf '%-18s %12s %14s %10s\n' SET INTERVALS BP MB
  for n in target bait; do
    awk -v n="$n" '{c++; s+=$3-$2} END{printf "%-18s %12d %14d %10.2f\n", n, c, s, s/1e6}' \
      "$BED_DIR/${n}.norm.b37.bed"
  done
  echo
  echo "contigs: $(cut -f1 "$BED_DIR/target.norm.b37.bed" | sort -u | tr '\n' ' ')"
} | tee "$BED_DIR/capture_design.txt"

echo
echo ">>> depth denominator = the 'target' row above (expect ~35 Mb for SureSelect V7)"

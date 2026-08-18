#!/bin/bash
#SBATCH --job-name=fetch_b37
#SBATCH --output=fetch_somatic_b37_%j.out
#SBATCH --error=fetch_somatic_b37_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH -p compute
#SBATCH -A o250038
# =============================================================================
# Stage the Broad somatic-b37 resource bundle for the hg19 profile.
#
#   sbatch scripts/fetch_somatic_b37.sh          # or run directly on a node
#
# WHY THIS EXISTS. /common/db/human_ref/hg19/v0 is a Broad b37 bundle, but it
# carries only the *germline* resources (dbSNP, Mills, HapMap, Omni). None of
# the three Mutect2 somatic resources exist anywhere on this cluster in b37:
# af-only gnomAD, small_exac_common_3, and a Mutect2 panel of normals. Neither
# does any b37 ClairS-TO panel of normals. This script fetches the first three
# from the Broad public bucket and derives the ClairS-TO set from them.
#
# The bucket files are stamped
#     ##contig=<ID=1,length=249250621,assembly=Homo_sapiens_assembly19.fasta>
# i.e. the exact reference our BAMs are aligned to - no liftover, no renaming.
#
# NOTE: raptor has no bgzip and no tabix binary. bcftools does both jobs:
#     bcftools view -Oz   ==  bgzip
#     bcftools index -t   ==  tabix -p vcf
# Do not "simplify" this back to bgzip/tabix; it will not run here.
#
# Idempotent: anything already present with a matching index is skipped, so a
# failed run can simply be resubmitted.
# =============================================================================
set -euo pipefail

BUNDLE_URL=${BUNDLE_URL:-https://storage.googleapis.com/gatk-best-practices/somatic-b37}
OUT_DIR=${OUT_DIR:-/common/db/human_ref/hg19/somatic_b37}
CLAIRS_DIR="$OUT_DIR/clairs_to"
REF_FAI=${REF_FAI:-/common/db/human_ref/hg19/v0/Homo_sapiens_assembly19.fasta.fai}
DBSNP_B37=${DBSNP_B37:-/common/db/human_ref/hg19/v0/dbsnp_138.b37.vcf.gz}
THREADS=${SLURM_CPUS_PER_TASK:-4}

module load bcftools/1.21 2>/dev/null || true
command -v bcftools >/dev/null || { echo "bcftools not on PATH" >&2; exit 1; }

mkdir -p "$CLAIRS_DIR"

# ---------------------------------------------------------------------------
# fetch <remote-plain-vcf> <local-bgzipped-name>
# Downloads to .part, verifies the byte count against Content-Length, and only
# then compresses + indexes. A truncated 14 GB gnomAD must never be promoted to
# a usable-looking resource.
# ---------------------------------------------------------------------------
fetch() {
  local remote=$1 out=$2
  if [[ -s "$OUT_DIR/$out" && -s "$OUT_DIR/$out.tbi" ]]; then
    echo "[skip ] $out (already staged)"
    return
  fi

  local part="$OUT_DIR/.part.$remote"
  local expect got
  expect=$(curl -sSIL "$BUNDLE_URL/$remote" \
           | awk 'BEGIN{IGNORECASE=1}/^content-length:/{v=$2}END{gsub(/\r/,"",v);print v}')
  echo "[get  ] $remote (expecting ${expect:-unknown} bytes)"
  curl -fSL --retry 5 --retry-delay 10 -o "$part" "$BUNDLE_URL/$remote"

  got=$(stat -c%s "$part")
  if [[ -n "$expect" && "$got" != "$expect" ]]; then
    echo "ERROR: $remote truncated - got $got bytes, expected $expect" >&2
    rm -f "$part"
    exit 1
  fi

  echo "[bgzip] $out"
  bcftools view --threads "$THREADS" -Oz -o "$OUT_DIR/.part.$out" "$part"
  bcftools index --threads "$THREADS" -t -f "$OUT_DIR/.part.$out"
  mv -f "$OUT_DIR/.part.$out"     "$OUT_DIR/$out"
  mv -f "$OUT_DIR/.part.$out.tbi" "$OUT_DIR/$out.tbi"
  rm -f "$part"
}

fetch af-only-gnomad.raw.sites.vcf  af-only-gnomad.raw.sites.b37.vcf.gz
fetch small_exac_common_3.vcf       small_exac_common_3.b37.vcf.gz
fetch Mutect2-exome-panel.vcf       Mutect2-exome-panel.b37.vcf.gz

# ---------------------------------------------------------------------------
# ClairS-TO panel of normals, b37.
#
# The four PoNs shipped with ClairS-TO (/project/o240001_SBUFF67/ClairS-TO/PON)
# are GRCh38 and carry no index at all. Three of them have a direct b37
# counterpart, built here as sites-only VCFs. The fourth, CoLoRSdb, is
# GRCh38-only: its coordinates would need a genuine liftover, not a contig
# rename, so it is deliberately omitted. main.nf accepts a 3-entry PoN list.
# ---------------------------------------------------------------------------
derive() { # <out-name> <in-vcf> [extra bcftools view args...]
  local out=$1; shift
  local in=$1;  shift
  if [[ -s "$CLAIRS_DIR/$out" && -s "$CLAIRS_DIR/$out.tbi" ]]; then
    echo "[skip ] clairs_to/$out (already built)"
    return
  fi
  echo "[build] clairs_to/$out"
  bcftools view --threads "$THREADS" -G "$@" -Oz -o "$CLAIRS_DIR/.part.$out" "$in"
  bcftools index --threads "$THREADS" -t -f "$CLAIRS_DIR/.part.$out"
  mv -f "$CLAIRS_DIR/.part.$out"     "$CLAIRS_DIR/$out"
  mv -f "$CLAIRS_DIR/.part.$out.tbi" "$CLAIRS_DIR/$out.tbi"
}

derive gnomad.af-ge-0.001.b37.sites.vcf.gz \
       "$OUT_DIR/af-only-gnomad.raw.sites.b37.vcf.gz" -i 'INFO/AF>=0.001'
derive dbsnp.b138.b37.sites.vcf.gz  "$DBSNP_B37"
derive 1000g-pon.b37.sites.vcf.gz   "$OUT_DIR/Mutect2-exome-panel.b37.vcf.gz"

# ---------------------------------------------------------------------------
# Normalise contig headers.
#
# af-only-gnomad.raw.sites.vcf ships from the Broad bucket with NO ##contig
# header lines at all. bcftools synthesises bare `##contig=<ID=1>` entries for
# whichever contigs it happens to encounter - 23 of them, no length=, no Y, no
# MT - and that thin dictionary propagates into everything derived from it.
# Rewrite the header from the reference .fai so every file carries the full
# 85-contig dictionary with lengths, matching the reference exactly.
#
# reheader only rewrites the header block, so this is cheap; it is also
# idempotent, guarded on whether length= is already present.
# ---------------------------------------------------------------------------
normalize_header() {
  local f=$1
  if [[ $(bcftools view -h "$f" | grep -c '^##contig=<ID=[^,>]*,length=') -gt 0 ]]; then
    echo "[skip ] header already has contig lengths: $(basename "$f")"
    return
  fi
  echo "[rehdr] $(basename "$f")"
  bcftools reheader --fai "$REF_FAI" -o "$f.rehdr" "$f"
  mv -f "$f.rehdr" "$f"
  bcftools index --threads "$THREADS" -t -f "$f"
}

for f in "$OUT_DIR"/*.vcf.gz "$CLAIRS_DIR"/*.vcf.gz; do
  normalize_header "$f"
done

# ---------------------------------------------------------------------------
# Contig assertion. Every contig carrying records must exist in the reference
# .fai. A mismatch here is the failure mode that produces silent empty output
# downstream (Mutect2 emitting nothing, GetPileupSummaries writing 0 rows).
# ---------------------------------------------------------------------------
REF_CONTIGS=$(mktemp)
trap 'rm -f "$REF_CONTIGS"' EXIT
cut -f1 "$REF_FAI" | sort -u > "$REF_CONTIGS"

check_contigs() {
  local f=$1 bad
  bad=$(bcftools index --stats "$f" | cut -f1 | sort -u | grep -vxFf "$REF_CONTIGS" || true)
  if [[ -n "$bad" ]]; then
    echo "ERROR: $f carries contigs absent from $REF_FAI:" >&2
    echo "$bad" >&2
    exit 1
  fi
  printf '[ok   ] %-46s %s records\n' "$(basename "$f")" \
    "$(bcftools index --stats "$f" | awk -F'\t' '{s+=$3} END{print s+0}')"
}

echo
echo "=== contig check against $REF_FAI ==="
for f in "$OUT_DIR"/*.vcf.gz "$CLAIRS_DIR"/*.vcf.gz; do
  check_contigs "$f"
done

echo
echo "b37 somatic bundle staged under $OUT_DIR"

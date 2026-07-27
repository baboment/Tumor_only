#!/bin/bash

set -euo pipefail

WORKDIR="${1:-/project/o240001_SBUFF67/script_zam/Tumor_only/work/e4/16007d01bd978de7357f0240c23aa3}"

if [[ ! -d "$WORKDIR" ]]; then
  echo "Work dir not found: $WORKDIR" >&2
  exit 1
fi

echo "Work dir: $WORKDIR"
echo

for f in .command.begin .command.run .command.sh .command.log .command.err .exitcode; do
  if [[ -f "${WORKDIR}/${f}" ]]; then
    echo "=== ${f} ==="
    ls -lh "${WORKDIR}/${f}"
    echo
  fi
done

if [[ -f "${WORKDIR}/.command.log" ]]; then
  echo "=== tail .command.log ==="
  tail -n 50 "${WORKDIR}/.command.log"
  echo
fi

if [[ -f "${WORKDIR}/.command.err" ]]; then
  echo "=== tail .command.err ==="
  tail -n 50 "${WORKDIR}/.command.err"
  echo
fi

if [[ -d "${WORKDIR}/clairs_to_out" ]]; then
  echo "=== clairs_to_out disk usage ==="
  du -sh "${WORKDIR}/clairs_to_out"
  echo

  echo "=== newest files under clairs_to_out ==="
  find "${WORKDIR}/clairs_to_out" -type f -printf '%TY-%Tm-%Td %TH:%TM:%TS %12s %p\n' | sort | tail -n 30
  echo
fi

echo "=== final outputs in work dir ==="
find "$WORKDIR" -maxdepth 2 -type f \( -name 'snv.vcf.gz' -o -name 'indel.vcf.gz' -o -name '*.clairs.vcf.gz' -o -name '*.tbi' \) -printf '%TY-%Tm-%Td %TH:%TM:%TS %12s %p\n' | sort || true

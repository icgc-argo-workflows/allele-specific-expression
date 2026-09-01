#!/usr/bin/env bash

set -euo pipefail

repository_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
smoke_root="$(mktemp -d "${TMPDIR:-/tmp}/pipeline-ase-smoke.XXXXXX")"
trap 'rm -rf "${smoke_root}"' EXIT

cd "${repository_root}"
nextflow run pipeline-ase.nf \
    --bam tests/data/smoke/smoke.sam \
    --vcf tests/data/smoke/smoke.vcf \
    --fa tests/data/smoke/smoke.fa \
    --out_dir "${smoke_root}/out" \
    --min_depth 1 \
    -work-dir "${smoke_root}/work" \
    -ansi-log false

test -s "${smoke_root}/out/smoke.tsv"
test -s "${smoke_root}/out/smoke.ase.log"
test -s "${smoke_root}/out/smoke.vaf.png"
grep -q $'chr1\t100' "${smoke_root}/out/smoke.tsv"

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
    --gtf tests/data/smoke/smoke.gtf \
    --assembly smoke \
    --mapp tests/data/smoke/smoke.bedgraph \
    --out_dir "${smoke_root}/out" \
    --debug_out "${smoke_root}/debug" \
    --min_depth 1 \
    -work-dir "${smoke_root}/work" \
    -ansi-log false

test -s "${smoke_root}/out/smoke.tsv"
test -s "${smoke_root}/out/smoke.ase.log"
test -s "${smoke_root}/out/smoke.vaf.png"
test -s "${smoke_root}/out/smoke.gene_annotation.log"
test -s "${smoke_root}/out/smoke.hap.tsv"
test -s "${smoke_root}/out/smoke.hap.log"
test -s "${smoke_root}/out/software_versions.yml"
test -s "${smoke_root}/out/execution_trace.txt"
test -s "${smoke_root}/out/execution_report.html"
test -s "${smoke_root}/out/execution_timeline.html"
test -s "${smoke_root}/debug/smoke.filtered.vcf"
grep -q $'chr1\t100' "${smoke_root}/out/smoke.tsv"
grep -q 'SMOKE_GENE' "${smoke_root}/out/smoke.tsv"
grep -q $'SMOKE_GENE\t100' "${smoke_root}/out/smoke.hap.tsv"
if grep -q $'chr1\t120' "${smoke_root}/debug/smoke.filtered.vcf"; then
    echo "Biallelic indel was not removed from the filtered VCF." >&2
    exit 1
fi

#!/usr/bin/env bash

set -euo pipefail

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
reference_directory="${script_directory}/refs"
mkdir -p "${reference_directory}"

verify_file() {
    local destination="$1"
    local expected_sha256="${2:-}"
    local actual_sha256
    actual_sha256="$(sha256sum "${destination}" | cut -d ' ' -f 1)"
    if [[ -n "${expected_sha256}" && "${actual_sha256}" != "${expected_sha256}" ]]; then
        echo "Checksum mismatch for ${destination}." >&2
        return 1
    fi

    local checksum_file="${destination}.sha256"
    if [[ -s "${checksum_file}" ]]; then
        local recorded_sha256
        recorded_sha256="$(cut -d ' ' -f 1 "${checksum_file}")"
        if [[ "${actual_sha256}" != "${recorded_sha256}" ]]; then
            echo "Downloaded file changed since its checksum was recorded: ${destination}." >&2
            return 1
        fi
    else
        printf '%s  %s\n' "${actual_sha256}" "$(basename "${destination}")" > "${checksum_file}"
    fi
}

download_file() {
    local url="$1"
    local destination="$2"
    local expected_sha256="${3:-}"

    if [[ ! -s "${destination}" ]]; then
        local partial_file="${destination}.part"
        curl \
            --fail \
            --location \
            --retry 5 \
            --retry-all-errors \
            --continue-at - \
            --output "${partial_file}" \
            "${url}"
        mv "${partial_file}" "${destination}"
    fi

    verify_file "${destination}" "${expected_sha256}"
}

reference_fasta="${reference_directory}/GRCh38_Verily_v1.genome.fa"
reference_bgzip="${reference_fasta}.gz"
if [[ ! -s "${reference_bgzip}" ]]; then
    download_file \
        "https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa" \
        "${reference_fasta}" \
        "${GRCH38_SHA256:-}"
    bgzip --force "${reference_fasta}"
fi
bgzip --test "${reference_bgzip}"
verify_file "${reference_bgzip}" "${GRCH38_BGZIP_SHA256:-}"

download_file \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz" \
    "${reference_directory}/gencode.v41.annotation.gtf.gz" \
    "${GENCODE_V41_SHA256:-}"
gzip --test "${reference_directory}/gencode.v41.annotation.gtf.gz"

download_file \
    "https://bismap.hoffmanlab.org/raw/hg38/k50.umap.bedgraph.gz" \
    "${reference_directory}/k50.umap.bedgraph.gz" \
    "${UMAP_K50_SHA256:-}"
gzip --test "${reference_directory}/k50.umap.bedgraph.gz"

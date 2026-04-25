#!/usr/bin/env bash
# Download the 5 Tarnita 2023 unannotated genome assemblies from NCBI.
# Run from: /home/pmt5gt/data/CarnivorEnzyme
# Requires: ncbi-datasets-cli (conda install -c conda-forge ncbi-datasets-cli)

set -euo pipefail

GENOME_DIR="resources/genomes/tarnita2023"
mkdir -p "$GENOME_DIR"

declare -A GENOMES=(
    ["darlingtonia_californica"]="GCA_032270445.1"
    ["heliamphora_ciliata"]="GCA_032360265.1"
    ["pinguicula_moranensis"]="GCA_028565015.1"
    ["drosera_spatulata"]="GCA_035668315.1"
    ["sarracenia_alata"]="GCA_019775975.1"
)

for species in "${!GENOMES[@]}"; do
    acc="${GENOMES[$species]}"
    out_fa="${GENOME_DIR}/${species}.fna"

    if [[ -f "$out_fa" ]]; then
        echo "[SKIP] $species already downloaded"
        continue
    fi

    echo "[DOWNLOAD] $species ($acc)"
    datasets download genome accession "$acc" \
        --include genome \
        --filename "${GENOME_DIR}/${species}_tmp.zip"

    unzip -p "${GENOME_DIR}/${species}_tmp.zip" \
        "ncbi_dataset/data/${acc}/*.fna" \
        > "$out_fa"

    rm -f "${GENOME_DIR}/${species}_tmp.zip"

    n_contigs=$(grep -c '^>' "$out_fa")
    size_mb=$(du -sm "$out_fa" | cut -f1)
    echo "[OK] $species: ${n_contigs} contigs, ${size_mb} MB"
done

echo "All genomes downloaded to $GENOME_DIR"

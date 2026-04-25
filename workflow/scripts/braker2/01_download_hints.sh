#!/usr/bin/env bash
# Download protein hint sequences for BRAKER2.
# Uses Arabidopsis thaliana (TAIR10) + Cephalotus follicularis (Fukushima 2017)
# as the evidence proteomes — both are well-annotated and closely enough related
# to produce reliable gene model evidence for carnivorous plant genomes.
# Run from: /home/pmt5gt/data/CarnivorEnzyme

set -euo pipefail

HINTS_DIR="resources/braker2_hints"
mkdir -p "$HINTS_DIR"

# --- Arabidopsis thaliana TAIR10 proteome ---
ATHAL_FA="${HINTS_DIR}/arabidopsis_thaliana.faa"
if [[ ! -f "$ATHAL_FA" ]]; then
    echo "[DOWNLOAD] Arabidopsis thaliana proteome (TAIR10)"
    datasets download genome accession GCF_000001735.4 \
        --include protein \
        --filename "${HINTS_DIR}/athal_tmp.zip"
    unzip -p "${HINTS_DIR}/athal_tmp.zip" \
        "ncbi_dataset/data/GCF_000001735.4/protein.faa" \
        > "$ATHAL_FA"
    rm -f "${HINTS_DIR}/athal_tmp.zip"
    n=$(grep -c '^>' "$ATHAL_FA")
    echo "[OK] Arabidopsis: $n proteins"
else
    echo "[SKIP] Arabidopsis already downloaded"
fi

# --- Cephalotus follicularis (Fukushima 2017 genome, GCA_001941015.1) ---
CEPH_FA="${HINTS_DIR}/cephalotus_follicularis.faa"
if [[ ! -f "$CEPH_FA" ]]; then
    echo "[DOWNLOAD] Cephalotus follicularis proteome"
    datasets download genome accession GCA_001941015.1 \
        --include protein \
        --filename "${HINTS_DIR}/ceph_tmp.zip"
    unzip -p "${HINTS_DIR}/ceph_tmp.zip" \
        "ncbi_dataset/data/GCA_001941015.1/protein.faa" \
        > "$CEPH_FA"
    rm -f "${HINTS_DIR}/ceph_tmp.zip"
    n=$(grep -c '^>' "$CEPH_FA")
    echo "[OK] Cephalotus: $n proteins"
else
    echo "[SKIP] Cephalotus already downloaded"
fi

# --- Dionaea muscipula (Palfalvi 2020, GCA_032404045.1) ---
# Includes nearest relative of Darlingtonia/Heliamphora/Drosera/Pinguicula
DION_FA="${HINTS_DIR}/dionaea_muscipula.faa"
if [[ ! -f "$DION_FA" ]]; then
    echo "[DOWNLOAD] Dionaea muscipula proteome"
    datasets download genome accession GCA_032404045.1 \
        --include protein \
        --filename "${HINTS_DIR}/dion_tmp.zip"
    unzip -p "${HINTS_DIR}/dion_tmp.zip" \
        "ncbi_dataset/data/GCA_032404045.1/protein.faa" \
        > "$DION_FA"
    rm -f "${HINTS_DIR}/dion_tmp.zip"
    n=$(grep -c '^>' "$DION_FA")
    echo "[OK] Dionaea: $n proteins"
else
    echo "[SKIP] Dionaea already downloaded"
fi

# --- Combine into single hints file ---
COMBINED="${HINTS_DIR}/combined_hints.faa"
cat "$ATHAL_FA" "$CEPH_FA" "$DION_FA" > "$COMBINED"
total=$(grep -c '^>' "$COMBINED")
echo "[OK] Combined hints: $total proteins → $COMBINED"

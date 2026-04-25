#!/usr/bin/env bash
# Submit one BRAKER2 SLURM job per Tarnita 2023 species.
# Run from: /home/pmt5gt/data/CarnivorEnzyme
# Prerequisites: 00_download_genomes.sh and 01_download_hints.sh must have completed.

set -euo pipefail

GENOME_DIR="resources/genomes/tarnita2023"
HINTS_FA="resources/braker2_hints/combined_hints.faa"
LOG_DIR="logs/braker2"
TEMPLATE="workflow/scripts/braker2/braker2_job.slurm"

mkdir -p "$LOG_DIR"

# Verify prerequisites
if [[ ! -f "$HINTS_FA" ]]; then
    echo "ERROR: hints file not found: $HINTS_FA"
    echo "Run 01_download_hints.sh first"
    exit 1
fi

declare -A GENOMES=(
    ["darlingtonia_californica"]="GCA_032270445.1"
    ["heliamphora_ciliata"]="GCA_032360265.1"
    ["pinguicula_moranensis"]="GCA_028565015.1"
    ["drosera_spatulata"]="GCA_035668315.1"
    ["sarracenia_alata"]="GCA_019775975.1"
)

SUBMITTED=()

for species in "${!GENOMES[@]}"; do
    genome_fa="${GENOME_DIR}/${species}.fna"

    if [[ ! -f "$genome_fa" ]]; then
        echo "[SKIP] $species — genome not found: $genome_fa"
        continue
    fi

    # Skip if BRAKER2 output already exists
    if [[ -f "results/braker2/${species}/augustus.hints.aa" ]]; then
        echo "[SKIP] $species — already has protein output"
        continue
    fi

    # Fill in the template
    job_script="${LOG_DIR}/braker2_${species}.slurm"
    sed \
        -e "s|SPECIES_PLACEHOLDER|${species}|g" \
        -e "s|GENOME_FA_PLACEHOLDER|${genome_fa}|g" \
        "$TEMPLATE" > "$job_script"

    job_id=$(sbatch --parsable "$job_script")
    echo "[SUBMITTED] $species → job $job_id"
    SUBMITTED+=("$species:$job_id")
done

echo ""
echo "=== Submitted ${#SUBMITTED[@]} jobs ==="
for entry in "${SUBMITTED[@]}"; do
    echo "  $entry"
done
echo ""
echo "Monitor with: squeue -u pmt5gt"
echo "After completion, run: bash workflow/scripts/braker2/03_extract_proteins.sh"

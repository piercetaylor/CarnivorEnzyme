#!/usr/bin/env bash
#SBATCH --job-name=braker2_SPECIES
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/braker2/braker2_SPECIES_%j.out
#SBATCH --error=logs/braker2/braker2_SPECIES_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pmt5gt@virginia.edu

# BRAKER2 gene prediction with protein hints.
# This template is filled in by 02_submit_all.sh — do not run directly.
# SPECIES and GENOME_FA are substituted at submission time.

set -euo pipefail

SPECIES="SPECIES_PLACEHOLDER"
GENOME_FA="GENOME_FA_PLACEHOLDER"
HINTS_FA="resources/braker2_hints/combined_hints.faa"
WORKDIR="results/braker2/${SPECIES}"
AUGUSTUS_CONFIG_PATH="${HOME}/augustus_config"

module load miniconda3 2>/dev/null || true
conda activate carnivorenzyme

mkdir -p "$WORKDIR" logs/braker2

echo "=== BRAKER2: $SPECIES ==="
echo "Genome: $GENOME_FA"
echo "Hints:  $HINTS_FA"
echo "Output: $WORKDIR"
echo "Started: $(date)"

# Copy AUGUSTUS config to writable location if needed
if [[ ! -d "$AUGUSTUS_CONFIG_PATH" ]]; then
    cp -r "$(which augustus | xargs dirname)/../config" "$AUGUSTUS_CONFIG_PATH"
fi
export AUGUSTUS_CONFIG_PATH

braker.pl \
    --genome="$GENOME_FA" \
    --prot_seq="$HINTS_FA" \
    --softmasking \
    --cores=16 \
    --workingdir="$WORKDIR" \
    --gff3 \
    --prg=gth \
    --species="carnivor_${SPECIES}" \
    --AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH"

echo "=== BRAKER2 complete: $SPECIES ==="
echo "Finished: $(date)"

# Quick stats
if [[ -f "${WORKDIR}/braker.gtf" ]]; then
    n_genes=$(grep -c $'\tgene\t' "${WORKDIR}/braker.gtf" || echo 0)
    echo "Predicted genes: $n_genes"
fi

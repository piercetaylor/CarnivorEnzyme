#!/usr/bin/env bash
# After BRAKER2 jobs complete: collect per-species protein FASTAs and
# run length filtering to keep only full-length predictions (>50 aa).
# Run from: /home/pmt5gt/data/CarnivorEnzyme

set -euo pipefail

BRAKER_DIR="results/braker2"
OUT_DIR="resources/accessions/tarnita2023"

mkdir -p "$OUT_DIR"

SPECIES=(
    "darlingtonia_californica"
    "heliamphora_ciliata"
    "pinguicula_moranensis"
    "drosera_spatulata"
    "sarracenia_alata"
)

for species in "${SPECIES[@]}"; do
    braker_aa="${BRAKER_DIR}/${species}/augustus.hints.aa"

    if [[ ! -f "$braker_aa" ]]; then
        echo "[MISSING] $species — ${braker_aa} not found (job may still be running)"
        continue
    fi

    out_fa="${OUT_DIR}/${species}.faa"

    # Filter: keep proteins ≥50 aa (removes fragments from incomplete gene models)
    python3 - <<EOF
from Bio import SeqIO
import sys

records = list(SeqIO.parse("${braker_aa}", "fasta"))
kept = [r for r in records if len(r.seq) >= 50]
SeqIO.write(kept, "${out_fa}", "fasta")
print(f"${species}: {len(records)} predicted → {len(kept)} kept (≥50 aa)")
EOF

    if [[ $? -eq 0 ]]; then
        n=$(grep -c '^>' "$out_fa")
        echo "[OK] $species → $out_fa ($n sequences)"
    fi
done

echo ""
echo "Protein FASTAs collected in $OUT_DIR"
echo "Next step: run 04_run_hmmer_scan.sh to identify enzyme families"

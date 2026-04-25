#!/usr/bin/env bash
# Run hmmsearch against Pfam HMMs to identify target enzyme families
# in each BRAKER2-predicted proteome. Writes per-family hit FASTAs
# ready to be added to enzyme_families.yaml as local: entries.
#
# Run from: /home/pmt5gt/data/CarnivorEnzyme
# Requires: hmmer (conda install -c bioconda hmmer), pfam_scan or hmmfetch

set -euo pipefail

PROTEOME_DIR="resources/accessions/tarnita2023"
HMM_DIR="resources/hmm"
HIT_DIR="resources/accessions/tarnita2023_by_family"

mkdir -p "$HMM_DIR" "$HIT_DIR"

# --- Download Pfam HMMs if not present ---
PFAM_HMM="${HMM_DIR}/Pfam-A.hmm"
if [[ ! -f "$PFAM_HMM" ]]; then
    echo "[DOWNLOAD] Pfam-A.hmm (this may take ~10 min)"
    wget -q -O "${HMM_DIR}/Pfam-A.hmm.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    gunzip "${HMM_DIR}/Pfam-A.hmm.gz"
    hmmpress "$PFAM_HMM"
fi

# --- Define target HMM accessions per enzyme family ---
declare -A FAMILY_HMMS=(
    # GH19 chitinase — PF00182 (Glyco_hydro_19)
    ["chitinases_gh19"]="PF00182"
    # Purple acid phosphatase — PF00149 (Calcineurin-like) + PF16656 (PAP N-term)
    # Use PF00149 as the conserved catalytic domain marker
    ["purple_acid_phosphatase"]="PF00149"
    # RNase T2 — PF00445 (RNase T2 family)
    ["rnase_t2"]="PF00445"
    # A1B aspartic protease — PF00026 (Asp) + PF01951 (Asp_protease_2/PSI)
    ["nepenthesins"]="PF00026"
)

# Length cutoffs per family (minimum aa for a plausible full-length prediction)
declare -A MIN_LEN=(
    ["chitinases_gh19"]="200"
    ["purple_acid_phosphatase"]="350"
    ["rnase_t2"]="180"
    ["nepenthesins"]="280"
)

SPECIES=(
    "darlingtonia_californica"
    "heliamphora_ciliata"
    "pinguicula_moranensis"
    "drosera_spatulata"
    "sarracenia_alata"
)

for species in "${SPECIES[@]}"; do
    proteome="${PROTEOME_DIR}/${species}.faa"

    if [[ ! -f "$proteome" ]]; then
        echo "[SKIP] $species — proteome not found"
        continue
    fi

    echo "=== $species ==="

    for family in "${!FAMILY_HMMS[@]}"; do
        pfam_acc="${FAMILY_HMMS[$family]}"
        min_len="${MIN_LEN[$family]}"
        hit_tbl="${HIT_DIR}/${species}_${family}_hits.tbl"
        hit_fa="${HIT_DIR}/${species}_${family}.faa"

        # Extract individual HMM from Pfam-A
        single_hmm="${HMM_DIR}/${pfam_acc}.hmm"
        if [[ ! -f "$single_hmm" ]]; then
            hmmfetch "$PFAM_HMM" "$pfam_acc" > "$single_hmm"
        fi

        # Run hmmsearch
        hmmsearch \
            --tblout "$hit_tbl" \
            -E 1e-5 \
            --cpu 4 \
            "$single_hmm" \
            "$proteome" \
            > /dev/null

        # Extract hit sequence IDs (col 1, skip comment lines)
        hit_ids=$(grep -v '^#' "$hit_tbl" | awk '{print $1}' | sort -u)

        if [[ -z "$hit_ids" ]]; then
            echo "  $family: 0 hits"
            continue
        fi

        # Extract sequences and filter by minimum length
        python3 - <<EOF
from Bio import SeqIO

proteome = "${proteome}"
hit_ids = set("""${hit_ids}""".strip().split())
min_len = ${min_len}
out_fa = "${hit_fa}"

records = {r.id: r for r in SeqIO.parse(proteome, "fasta")}
hits = [records[hid] for hid in hit_ids if hid in records]
kept = [r for r in hits if len(r.seq) >= min_len]

# Rename IDs to include species prefix
for r in kept:
    r.id = f"${species}|{r.id}"
    r.description = f"family=${family} species=${species} source=braker2"

SeqIO.write(kept, out_fa, "fasta")
print(f"  ${family}: {len(hit_ids)} hmmer hits → {len(kept)} kept (≥{min_len} aa) → {out_fa}")
EOF

    done
done

echo ""
echo "=== Summary ==="
echo "Family hit FASTAs written to $HIT_DIR"
echo ""
echo "Next steps:"
echo "  1. Inspect hits manually (check lengths, review annotations)"
echo "  2. Add local: entries to enzyme_families.yaml:"
echo "     e.g.  Darlingtonia_californica:"
echo "             - local:resources/accessions/tarnita2023_by_family/darlingtonia_californica_chitinases_gh19.faa"
echo "  3. Re-run fetch_sequences.py (after adding local: support)"

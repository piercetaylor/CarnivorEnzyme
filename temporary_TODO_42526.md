# Manual TODOs — CarnivorEnzyme (April 2026)

These cannot be resolved by writing code — they require manual data retrieval or
HPC job submission. Sorted by priority (highest-impact first).

---

## 1. BRAKER2 — Annotate Tarnita 2023 Genomes (15 missing accessions)

**Why:** 5 species have raw genome assemblies on NCBI but no protein annotation.
Proteins are needed for GH19, PAP, RNase T2, and A1B AP enzyme families.

**Status:** SLURM scripts written → `workflow/scripts/braker2/`
**Location to run:** `/home/pmt5gt/data/CarnivorEnzyme` on Hellbender

| Species | Assembly accession | Missing families | Priority |
|---------|--------------------|-----------------|----------|
| *Darlingtonia californica* | GCA_032270445.1 | GH19, PAP, RNase T2, AP | HIGH |
| *Heliamphora ciliata* | GCA_032360265.1 | GH19, PAP, RNase T2, AP | HIGH |
| *Pinguicula moranensis* | GCA_028565015.1 | GH19, PAP, RNase T2, AP | HIGH |
| *Drosera spatulata* | GCA_035668315.1 | GH19, AP | MEDIUM |
| *Sarracenia alata* | GCA_019775975.1 | GH19 | LOW (S. purpurea covers Ericales) |

**Steps:**
```bash
# 1. On Hellbender login node:
cd /home/pmt5gt/data/CarnivorEnzyme
bash workflow/scripts/braker2/00_download_genomes.sh      # ~30 min
bash workflow/scripts/braker2/01_download_hints.sh        # ~20 min
bash workflow/scripts/braker2/02_submit_all.sh            # submits 5 SLURM jobs
# 3. Wait ~12-24 h for jobs to complete
bash workflow/scripts/braker2/03_extract_proteins.sh      # collect augustus.hints.aa
bash workflow/scripts/braker2/04_run_hmmer_scan.sh        # identify enzyme families in output
```

**After completion:** Add discovered accessions (as `local:` paths) to `enzyme_families.yaml`
and re-run `fetch_sequences.py`.

---

## 2. ENA Download — *Nepenthes mirabilis* Proteins (3 missing accessions)

**Why:** Goh 2025 (PLOS One, PRJEB86749) deposited the N. mirabilis genome at ENA only —
not submitted to NCBI Protein. Gene IDs already identified from S3 InterProScan file.

**Target gene IDs:**

| Family | Gene ID | Length | Domain evidence |
|--------|---------|--------|----------------|
| GH19 | Nep_mir_jg4646.t1 | 276 aa | PF00182 E=5.2e-43 |
| GH19 | Nep_mir_jg6978.t1 | 360 aa | Chitinase PANTHER |
| Neprosin | Nep_mir_jg11529.t1 | ~400 aa (full len TBD) | PF03080 E=7.8e-10 |
| Neprosin | Nep_mir_jg39601.t1 | ~400 aa (full len TBD) | DUF239 PANTHER |
| PAP | Nep_mir_jg12124.t1 | 426 aa | PF16656 + cd00839 + PANTHER PTHR22953 |
| PAP | Nep_mir_jg35245.t1 | 615 aa | cd00839 + PF14008 |
| PAP | Nep_mir_jg38588.t1 | 661 aa | PF16656 + PF14008 |

**Steps:**
```bash
# Download full protein FASTA from ENA
curl -L "https://www.ebi.ac.uk/ena/browser/api/fasta/PRJEB86749?download=true" \
    -o resources/accessions/n_mirabilis_all_proteins.fasta

# Index it
samtools faidx resources/accessions/n_mirabilis_all_proteins.fasta

# Extract targets
samtools faidx resources/accessions/n_mirabilis_all_proteins.fasta \
    Nep_mir_jg4646.t1 Nep_mir_jg6978.t1 \
    > resources/accessions/n_mirabilis_gh19.fasta

samtools faidx resources/accessions/n_mirabilis_all_proteins.fasta \
    Nep_mir_jg11529.t1 Nep_mir_jg39601.t1 \
    > resources/accessions/n_mirabilis_neprosin.fasta

samtools faidx resources/accessions/n_mirabilis_all_proteins.fasta \
    Nep_mir_jg12124.t1 Nep_mir_jg35245.t1 Nep_mir_jg38588.t1 \
    > resources/accessions/n_mirabilis_pap.fasta
```

**After completion:** Update `fetch_sequences.py` to accept `local:<path>` prefix, then
update YAML entries from `TODO` to `local:resources/accessions/n_mirabilis_*.fasta`.

---

## 3. JGI Download — *Utricularia gibba* Proteins (4 missing accessions)

**Why:** Only plastid proteins exist for U. gibba in NCBI/UniProt (confirmed).
Lan 2017 gene models are at JGI (GCA_002189035.1) in unitig format.

**Known gene IDs:**
- RNase T2: `unitig_26.g10301.t1` (Plant Methods 2020, PMC7149871)
- GH19: 5-copy class IV expansion confirmed in Lan 2017 — unitig IDs unknown, need HMM scan

**Steps:**
```bash
# Download Lan 2017 U. gibba protein FASTA from JGI
# URL: https://genome.jgi.doe.gov/portal/UtrgibA2/Utrgib1_proteins.aa.fasta.gz
# Requires JGI login. Download manually via browser or:
curl -b cookies.txt \
    "https://genome.jgi.doe.gov/portal/Utrgib1/download/Utrgib1_proteins.aa.fasta.gz" \
    -o resources/accessions/u_gibba_jgi_proteins.fasta.gz
gunzip resources/accessions/u_gibba_jgi_proteins.fasta.gz

# Run HMM scans for target families
hmmbuild resources/hmm/gh19.hmm resources/hmm/PF00182.stockholm
hmmbuild resources/hmm/metallophos.hmm resources/hmm/PF00149.stockholm
hmmbuild resources/hmm/rnase_t2.hmm resources/hmm/cd08148.stockholm
hmmbuild resources/hmm/asp_protease.hmm resources/hmm/PF00026.stockholm

for family in gh19 metallophos rnase_t2 asp_protease; do
    hmmsearch --tblout resources/accessions/u_gibba_${family}_hits.txt \
        -E 1e-5 resources/hmm/${family}.hmm \
        resources/accessions/u_gibba_jgi_proteins.fasta
done
```

**After completion:** Extract top hits per family; store as `resources/accessions/u_gibba_*.fasta`.

---

## 4. Permanent Gaps — Manuscript Methods Note Required

These cannot be resolved (no public data exists). Document in manuscript Methods:

| Gap | Details |
|-----|---------|
| *Byblis filifolia* GH19 | No genome published as of April 2026 |
| *N. rafflesiana* neprosin | No genome; *N. × ventrata* (C0HLV2) is hybrid surrogate |
| *N. bicalcarata* neprosin | No genome; mutualistic species, may have reduced neprosin |
| *S. purpurea* RNase T2 | Transcriptome (PRJNA80051) exists but no deposited accession |
| Aldrovanda GH19, PAP, A1B AP | Palfalvi 2020 nuclear genome NOT submitted to NCBI Protein |

**Suggested wording for Methods:**
> "Sequences for [X species] were unavailable in public repositories as of April 2026 due to [reason].
> These species are omitted from the [family] alignment and do not affect the convergence analysis,
> which requires a minimum of three independent origins represented."

---

## 5. YAML Updates After Manual Steps Complete

Once items 1–3 above are done, update `enzyme_families.yaml`:

- [ ] Replace `TODO` entries for 5 Tarnita 2023 species with `local:` paths to BRAKER2 output FASTAs
- [ ] Replace `TODO` entries for *N. mirabilis* GH19/neprosin/PAP with `local:` paths
- [ ] Replace `TODO` entries for *U. gibba* with `local:` paths
- [ ] Add `local:` support to `fetch_sequences.py` (copy files instead of fetching from API)
- [ ] Re-run `fetch_sequences.py` and verify new sequence counts

---

## 6. CDD Domain Verification (Low Priority — Pre-FoldX)

Before running FoldX, verify domain assignments for entries marked with CAUTION notes:

| Accession | Species | Concern |
|-----------|---------|---------|
| GMH02520.1, GMH01779.1 | *N. gracilis* PAP | Annotated "hypothetical protein" — verify Metallophos via CDD |
| GMG98952.1, GMG98423.1 | *N. gracilis* RNase T2 | Annotated "hypothetical protein" — verify RNase_T2_euk via CDD |
| KAL9268294.1, KAL9267125.1 | *D. capensis* RNase T2 | Verify cd08148 (RNase T2) domain |
| GAB2239788.1 | *D. rotundifolia* AP | Verify A1B subfamily (PLN03146) |
| BAW35441.1 | *S. purpurea* AP | Verify PSI insert presence for A1B placement |

**Tool:** `https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi` (web) or
`rpsblast -db /path/to/CDD -query accession.fasta -outfmt 6` (command line)

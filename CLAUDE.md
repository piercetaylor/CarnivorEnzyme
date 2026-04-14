# CLAUDE.md — CarnivorEnzyme Project Guide

> This document is the single source of truth for the CarnivorEnzyme project.
> It is designed to be placed at the root of the repository and read by Claude Code.
> Every method choice is justified against published literature.
> Every tool version, PDB code, and accession has been verified.

---

## 0. Project Summary

**Title:** CarnivorEnzyme — Structural Atlas of Convergently Evolved Digestive Enzymes in Carnivorous Plants

**Goal:** Predict 3D structures of digestive enzymes from all independently-evolved carnivorous plant lineages, map published convergent amino acid substitutions onto those structures, quantify thermodynamic consequences with FoldX, and predict substrate specificity shifts via molecular docking. Deliver a reproducible Snakemake pipeline, a searchable SQLite database, a Streamlit web interface, and a preprint-ready manuscript.

**Why this matters:** Convergent amino acid substitutions in carnivorous plant digestive enzymes are documented at the sequence level (Fukushima et al. 2017; Wakatake & Fukushima 2025) but their structural and functional consequences are unknown. The Fukushima lab's 2025 review explicitly flags this as an open question (Albert, Fukushima et al. 2025). No cross-lineage structural comparison of ANY carnivorous plant digestive enzyme family has been published.

---

## 1. Method Justification Table

Each tool was selected based on published benchmarks and practical constraints. Do not substitute tools without updating this table and re-verifying accuracy claims.

| Step | Tool | Version | Why This Tool | Key Benchmark / Citation | Known Limitations |
|------|------|---------|---------------|--------------------------|-------------------|
| Orthology | OrthoFinder | 3.x | Highest F-score on OrthoBench (8–33% better than OrthoMCL/OMA). Provides rooted gene trees and duplication events. | Emms & Kelly 2019, *Genome Biology* 20:238; Emms & Kelly 2015, *Genome Biology* 16:157 | Recall can drop for very distant species; mitigated by our lineages being within angiosperms |
| Alignment | MAFFT L-INS-i | 7.525+ | Most accurate mode for <200 sequences per family. Iterative refinement with local pairwise alignment. | Katoh & Standley 2013, *Mol. Biol. Evol.* 30:772 | Slow for >500 seqs; irrelevant here (expect <100 per family) |
| Alignment trimming | trimAl | 1.4.1 | Automated gap removal with `-automated1` flag. Removes columns with >50% gaps while preserving informative positions. | Capella-Gutiérrez et al. 2009, *Bioinformatics* 25:1972 | Aggressive trimming can remove real insertions; inspect alignments manually for Tier 1 families |
| Phylogeny | IQ-TREE2 | 2.3+ | ModelFinder Plus (`-m MFP`) for automatic model selection. UFBoot2 for ultrafast bootstrap. Consistently top-performing in simulation benchmarks. | Minh et al. 2020, *Mol. Biol. Evol.* 37:1530; Kalyaanamoorthy et al. 2017, *Nature Methods* 14:587 | UFBoot2 values are NOT directly comparable to standard bootstrap; treat ≥95 as strong support |
| Convergence detection | Ancestral reconstruction + counting | Custom script | Bayesian ancestral sequence reconstruction on per-family gene trees, then counting excess convergent substitutions on carnivorous lineage branches vs. background expectation. This is the method from Fukushima et al. 2017 and Fukushima & Pollock 2023. | Fukushima et al. 2017, *Nat. Ecol. Evol.* 1:0059; Fukushima & Pollock 2023, *Nat. Ecol. Evol.* | Requires correctly rooted gene trees; depends on OrthoFinder rooting quality |
| Structure prediction | **Chai-1** (primary) | 0.6.1+ (PyPI) | Apache 2.0 license, pip-installable, weights freely downloadable. Protein–ligand PoseBusters: 77% RMSD <2 Å, comparable to AF3 (76%). Single-sequence mode outperforms ESMFold. | Chai Discovery 2024, bioRxiv 2024.10.10.615955 | Requires GPU with bfloat16 (A100/H100/L40S recommended). Not as extensively validated as AF3 on all complex types. |
| Structure prediction (optional) | AlphaFold3 | 3.0.1 | Gold standard. Available on some HPC clusters (OSC, UMich, NIH). Requires weight access request from Google DeepMind — NOT guaranteed. | Abramson et al. 2024, *Nature* 630:493 | Weights restricted (CC-BY-NC-SA 4.0 + separate Terms of Use). Requires A100+ GPU. Check if Hellbender has access before planning around this. |
| Structure comparison | ABCFold | 1.0+ (PyPI) | Unified wrapper for AF3, Boltz-1, and Chai-1 with standardized I/O and comparison visualization. | Elliott et al. 2025, *Bioinformatics Advances* 5:vbaf153 | Does not add accuracy; adds convenience for multi-predictor comparison |
| Structural alignment | Foldseek | 9+ | 3Di+AA structural alphabet enables all-vs-all structural search orders of magnitude faster than TM-align. | van Kempen et al. 2024, *Nature Biotechnology* 42:243 | Approximate; use TM-align for final pairwise TM-scores on top hits |
| Structural alignment (pairwise) | TM-align | latest | Gold standard for structural superposition and TM-score calculation. | Zhang & Skolnick 2005, *Nucleic Acids Research* 33:2302 | Only pairwise; not scalable for all-vs-all on 200 structures |
| Solvent accessibility | DSSP | 4.x | Standard method for secondary structure and accessible surface area from PDB coordinates. | Kabsch & Sander 1983, *Biopolymers* 22:2577; updated Joosten et al. 2011 | Requires well-formed PDB; predicted structures from Chai-1 are compatible |
| Stability prediction | FoldX | 5.0+ | Fast empirical force field (~seconds/mutation). PositionScan for systematic scanning of convergent positions. Validated on >1000 mutants. | Schymkowitz et al. 2005, *Nucleic Acids Research* 33:W382; Guerois et al. 2002, *J. Mol. Biol.* 320:369 | **Destabilizing mutations identified correctly ~69% of the time; stabilizing only ~29%** (Buß et al. 2018). Appropriate for classifying convergent substitutions as destabilizing/neutral/stabilizing, but NOT for quantitative engineering. Academic license required (free registration). |
| Interface affinity | PRODIGY | web/CLI | Predicts binding affinity from intermolecular contacts. Validated on protein-protein benchmark 5.0. | Xue et al. 2016, *Bioinformatics* 32:3676 | Only for protein-protein interfaces, not protein-substrate |
| Molecular docking | AutoDock Vina | 1.2.5+ | Open-source, scriptable, well-validated for small molecules and short peptides (<15 residues). | Eberhardt et al. 2021, *J. Chem. Inf. Model.* 61:3891; Trott & Olson 2010, *J. Comput. Chem.* 31:455 | **Not suitable for peptides >15 residues.** The gliadin 33-mer is too large for rigid-body docking; use 8-mer fragments instead. Docking scores are relative rankings, not absolute affinities. |
| Receptor preparation | ADFR Suite | 1.0 | Replaces deprecated MGLTools/AutoDockTools for receptor and ligand preparation. Maintained by Scripps. | Ravindranath et al. 2015, *PLOS Comput. Biol.* 11:e1004586 | Requires separate installation (not conda-installable); Python 3 compatible |
| Ligand preparation | Meeko | 0.5+ | Generates PDBQT from PDB/SDF. Modern replacement for prepare_ligand4.py. pip-installable. | Forli lab, Scripps (GitHub: forlilab/Meeko) | Limited to small molecules and short peptides |
| Electrostatics | APBS | 3.4+ | Poisson-Boltzmann electrostatic calculations. Standard for comparing charge surfaces at different pH values. | Baker et al. 2001, *Proc. Natl. Acad. Sci.* 98:10037 | pH-dependent protonation states must be set upstream (use PDB2PQR with PROPKA) |
| Protonation | PDB2PQR + PROPKA | 3.6+ | Assigns protonation states at user-specified pH. Required upstream of APBS. | Dolinsky et al. 2004, *Nucleic Acids Research* 32:W665 | PROPKA pKa predictions have ~1 pH unit error; adequate for pH 2.5 vs pH 5 comparison |
| Expression quantification | Salmon | 1.10+ | Alignment-free, fast, accurate quasi-mapping for gene-level TPM. We only need relative expression for ~50 target genes, not genome-wide novel transcript discovery. | Patro et al. 2017, *Nature Methods* 14:417 | Requires a transcriptome index per species; genome annotation must exist |
| Workflow | Snakemake | 8.x | Python-native, conda integration, SLURM profile for HPC. Standard in computational biology. | Mölder et al. 2021, *F1000Research* 10:33 | — |
| Web interface | Streamlit | 1.x | Python-native, deploys to Streamlit Cloud. py3Dmol for 3D structure viewing in browser. | — | Not suitable for high-traffic production; fine for a community resource |
| Database | SQLite | 3.x | Single-file, zero-config, adequate for ~200 structures × ~50 measurements. No server needed. | — | No concurrent write access; irrelevant for this use case |

---

## 2. Experimental Structures Available for Validation

These PDB entries exist and must be used for benchmarking predictions. Do not predict structures for these — download them.

| PDB | Description | Resolution | Organism | Citation |
|-----|-------------|-----------|----------|----------|
| **7ZVA** | Neprosin zymogen (native) | 1.80 Å | *N. × ventrata* | Del Amo-Maestro et al. 2022, *Nat. Commun.* 13:4446 |
| **7ZU8** | Neprosin zymogen + Lu-Xo4 | 2.05 Å | *N. × ventrata* | Del Amo-Maestro et al. 2022 |
| **7ZVB** | Neprosin mature form | 2.35 Å | *N. × ventrata* | Del Amo-Maestro et al. 2022 |
| **7ZVC** | Neprosin mature (2nd crystal form) | 1.85 Å | *N. × ventrata* | Del Amo-Maestro et al. 2022 |

**No experimental structure exists for:**
- Any nepenthesin (crystallized at 2.8 Å by Kadek et al. 2016, *Acta Cryst. F* 72:24, but coordinates never deposited)
- Any droserasin
- Any dionain
- Any carnivorous plant chitinase, RNase T2, or purple acid phosphatase

---

## 3. Target Enzyme Families

### Tier 1 — Manuscript core (convergent substitutions documented in Fukushima 2017)

1. **Chitinases (GH19)** — convergent substitutions confirmed (Fukushima et al. 2017, Fig. 3)
2. **Purple acid phosphatases (PAPs)** — convergent substitutions confirmed (Fukushima et al. 2017; Nishimura et al. 2014)
3. **RNase T2 (S-like RNases)** — convergent substitutions confirmed (Fukushima et al. 2017)
4. **Aspartic proteases (nepenthesins/droserasins, A1B)** — primary digestive proteases, well-characterized biochemically
5. **Glutamic peptidases (neprosins, G3)** — crystal structures available (7ZVA-C), biotech applications active

### Tier 2 — Extended atlas (convergent substitutions not yet characterized)

6. Cysteine proteases (dionain, droserain)
7. Thaumatin-like proteins (TLPs)
8. β-1,3-glucanases (GH17)
9. Lipid transfer proteins
10. Esterases/lipases

**Build Tier 1 first. Add Tier 2 only after Tier 1 pipeline is validated and producing correct results.**

---

## 4. Carnivorous Plant Lineages (Independent Origins)

| Lineage | Order | Representative Species | Genome Available | Trap Type |
|---------|-------|----------------------|-----------------|-----------|
| *Nepenthes* | Caryophyllales | *N. gracilis* | Yes (decaploid, Saul et al. 2023, *Nat. Plants*) | Pitfall pitcher |
| *Drosera* | Caryophyllales | *D. capensis*, *D. spatulata* | Draft (*D. capensis*; Butts et al. 2016 for proteases) | Sticky trap |
| *Dionaea* | Caryophyllales | *D. muscipula* | Yes (tetraploid, Palfalvi et al. 2020) | Snap trap |
| *Aldrovanda* | Caryophyllales | *A. vesiculosa* | Partial | Snap trap (aquatic) |
| *Drosophyllum* | Caryophyllales | *D. lusitanicum* | No genome; transcriptome exists | Sticky trap |
| *Cephalotus* | Oxalidales | *C. follicularis* | Yes (Fukushima et al. 2017) | Pitfall pitcher |
| *Sarracenia* | Ericales | *S. purpurea* | Transcriptome only | Pitfall pitcher |
| *Utricularia* | Lamiales | *U. gibba* | Yes (Lan et al. 2017) | Suction bladder |
| *Pinguicula* | Lamiales | *P. vulgaris* | Partial | Sticky trap |

**Non-carnivorous outgroups** (required for ancestral reconstruction):
- *Arabidopsis thaliana* (Brassicales)
- *Vitis vinifera* (Vitales)
- *Solanum lycopersicum* (Solanales)

---

## 5. Repository Structure

```
CarnivorEnzyme/
├── CLAUDE.md                          # THIS FILE — project guide for Claude Code
├── README.md                          # Public-facing description
├── LICENSE                            # MIT
├── pyproject.toml                     # Python package config (PEP 621)
├── environment.yml                    # Primary conda environment
│
├── config/
│   ├── config.yaml                    # Pipeline parameters (see Section 7)
│   ├── enzyme_families.yaml           # Family definitions + accession lists
│   ├── species.yaml                   # Taxonomy, pitcher pH, genome accessions
│   └── substrates.yaml               # Docking substrate definitions
│
├── Snakefile                          # Master workflow entry point
├── workflow/
│   ├── rules/
│   │   ├── retrieve.smk               # NCBI/UniProt sequence retrieval
│   │   ├── orthology.smk              # OrthoFinder clustering
│   │   ├── alignment.smk              # MAFFT + trimAl
│   │   ├── phylogeny.smk              # IQ-TREE2
│   │   ├── convergence.smk            # Ancestral reconstruction + detection
│   │   ├── predict_structure.smk      # Chai-1 batch prediction
│   │   ├── structural_align.smk       # Foldseek + TM-align
│   │   ├── foldx.smk                  # RepairPDB + PositionScan
│   │   ├── docking.smk                # Vina substrate docking
│   │   ├── electrostatics.smk         # PDB2PQR + APBS
│   │   ├── expression.smk             # Salmon quantification
│   │   └── integrate.smk              # Build atlas SQLite + figures
│   │
│   ├── scripts/                       # One script per task, Click CLI
│   │   ├── fetch_sequences.py         # Entrez + UniProt batch retrieval
│   │   ├── detect_convergence.py      # Ancestral reconstruction + counting
│   │   ├── map_convergence.py         # Project convergent sites onto 3D
│   │   ├── predict_chai1.py           # Chai-1 batch wrapper
│   │   ├── assess_structure.py        # pLDDT/PAE extraction + QC
│   │   ├── classify_positions.py      # Active site / surface / interface
│   │   ├── run_foldx_repair.py        # FoldX RepairPDB
│   │   ├── run_foldx_scan.py          # FoldX PositionScan on convergent sites
│   │   ├── parse_foldx.py             # ΔΔG matrix extraction
│   │   ├── prepare_docking.py         # ADFR receptor + Meeko ligand prep
│   │   ├── run_docking.py             # Vina batch execution
│   │   ├── parse_docking.py           # Aggregate binding affinities
│   │   ├── run_electrostatics.py      # PDB2PQR → APBS pipeline
│   │   ├── quantify_expression.py     # Salmon wrapper
│   │   ├── build_atlas.py             # Assemble SQLite from all results
│   │   └── generate_figures.py        # Publication figures (matplotlib)
│   │
│   └── envs/                          # Per-rule conda environments
│       ├── bioinfo.yaml               # MAFFT, IQ-TREE2, OrthoFinder, trimAl
│       ├── structure.yaml             # Chai-1, Foldseek, TM-align, DSSP
│       ├── foldx.yaml                 # FoldX (mounted license)
│       ├── docking.yaml               # Vina, Meeko, OpenBabel
│       └── expression.yaml            # Salmon
│
├── resources/                         # Input data (gitignored except configs)
│   ├── accessions/                    # TSV files of NCBI/UniProt accessions
│   ├── structures/
│   │   └── experimental/              # 7ZVA.pdb, 7ZVB.pdb, 7ZVC.pdb, 7ZU8.pdb
│   └── substrates/                    # Docking ligand files (.pdb, .sdf)
│
├── results/                           # All pipeline output (gitignored)
│
├── webapp/
│   ├── app.py                         # Streamlit entry point
│   ├── pages/
│   │   ├── family_browser.py
│   │   ├── structure_viewer.py
│   │   ├── foldx_results.py
│   │   └── docking_results.py
│   └── requirements.txt
│
├── manuscript/
│   ├── main.tex
│   ├── figures/
│   └── references.bib
│
└── tests/
    ├── test_fetch.py
    ├── test_convergence.py
    ├── test_foldx_parse.py
    ├── test_docking_parse.py
    └── conftest.py                    # Shared fixtures
```

---

## 6. Conda Environment

```yaml
# environment.yml
name: carnivorenzyme
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - snakemake=8.20
  # --- Bioinformatics ---
  - mafft=7.525
  - iqtree=2.3.6
  - orthofinder=2.5.5
  - trimal=1.4.1
  - diamond=2.1                        # OrthoFinder uses DIAMOND by default
  - blast=2.16                         # Fallback for OrthoFinder
  # --- Structure ---
  - foldseek=9.427613d
  - dssp=4.4.0
  - openbabel=3.1.1
  # --- Electrostatics ---
  - apbs=3.4.1
  - pdb2pqr=3.6
  # --- Expression ---
  - salmon=1.10.3
  # --- Python core ---
  - biopython=1.84
  - pandas=2.2
  - numpy=1.26
  - scipy=1.13
  - matplotlib=3.9
  - seaborn=0.13
  - scikit-learn=1.5
  - click
  - pyyaml
  - tqdm
  - requests
  - pytest
  # --- pip-only ---
  - pip
  - pip:
    - chai_lab>=0.6.1                  # Structure prediction (requires GPU)
    - meeko>=0.5                       # Docking ligand preparation
    - vina>=1.2.5                      # AutoDock Vina Python bindings
    - prodigy-prot                     # Interface binding affinity
    - streamlit>=1.30
    - py3Dmol                          # 3D viewer for Streamlit
    - ABCFold>=1.0                     # Optional: multi-predictor comparison
```

### External tools NOT in conda (manual installation required)

| Tool | How to install | Notes |
|------|---------------|-------|
| FoldX 5.0+ | Register at foldxsuite.crg.eu → download → set `$FOLDX_BINARY` env var | Free for academic use. Binary only. Mount into container or set PATH. |
| ADFR Suite 1.0 | Download from ccsb.scripps.edu/adfr → install per instructions | Needed for receptor preparation. Python 3 compatible. |
| TM-align | Download from zhanggroup.org/TM-align → compile `TMalign.cpp` | Single C++ file, no dependencies. |

---

## 7. Build Order and Test Gates

Every phase has a validation checkpoint. Do not proceed past a gate until it passes.

### Phase 1: Sequence Retrieval & Orthology

**Scripts to build:**
1. `workflow/scripts/fetch_sequences.py` — Uses Biopython `Entrez` for NCBI and `requests` for UniProt REST API. Reads accession TSVs from `config/`. Writes one FASTA per species to `results/sequences/`. Rate-limits NCBI to 3 requests/second with API key.
2. Snakemake rule `retrieve.smk` — Maps accession configs to FASTA outputs.
3. Snakemake rule `orthology.smk` — Runs `orthofinder -f results/sequences/ -t 8 -M msa`. Outputs to `results/orthologs/`.
4. Snakemake rule `alignment.smk` — Per ortholog group: `mafft --localpair --maxiterate 1000 input.fa > aligned.fa`, then `trimal -in aligned.fa -out trimmed.fa -automated1`.
5. Snakemake rule `phylogeny.smk` — Per family: `iqtree2 -s trimmed.fa -m MFP -bb 1000 -nt AUTO`.

**Test gate 1:** Manually inspect nepenthesin orthogroup. Must contain Nep I + Nep II from *N. gracilis* (AB114914, AB114915), droserasins from *D. capensis*, and the *Cephalotus* aspartic protease from Fukushima 2017 supplementary data. Tree topology must place Caryophyllales nepenthesins/droserasins as monophyletic, with *Cephalotus* (Oxalidales) outgroup. If not, the orthology step has failed.

### Phase 2: Convergence Detection

**Scripts to build:**
6. `workflow/scripts/detect_convergence.py` — Implements: (a) ancestral sequence reconstruction using maximum likelihood on the per-family IQ-TREE gene tree, (b) identification of positions where ≥2 independent carnivorous lineage branches converge on the same derived amino acid, (c) statistical test comparing observed convergent count to expectation under a null model (Poisson, parameterized by branch lengths and substitution rates). Uses Biopython for tree parsing and ete3 or custom code for ancestral reconstruction. Output: TSV with columns `family, alignment_position, ancestral_aa, derived_aa, lineages, posterior_prob, p_value`.
7. Snakemake rule `convergence.smk`.

**Test gate 2:** Run on GH19 chitinases. Compare detected convergent positions to Fukushima et al. 2017, Fig. 3a. Must recover ≥80% of the published positions. If recall is <80%, the ancestral reconstruction or tree rooting is incorrect.

### Phase 3: Structure Prediction

**Scripts to build:**
8. `workflow/scripts/predict_chai1.py` — Reads FASTA of all unique sequences from ortholog groups. For each sequence: writes a FASTA input, calls `chai_lab.chai1.run_inference()` with `num_trunk_recycles=3`, `num_diffn_timesteps=200`, `seed=42` (5 seeds for Tier 1 families, 1 for Tier 2). Outputs .cif per prediction. Extracts pLDDT and PAE from output.
9. `workflow/scripts/assess_structure.py` — Reads Chai-1 output, computes per-residue pLDDT, flags regions with pLDDT <50 as disordered, 50-70 as low confidence, >70 as confident. Writes QC summary TSV.
10. `workflow/scripts/classify_positions.py` — For each convergent position from Phase 2: determine if it falls in (a) active site (within 5Å of known catalytic residues), (b) substrate-binding pocket (CASTp or fpocket), (c) surface-exposed (DSSP relative ASA >25%), (d) buried (DSSP relative ASA <5%), (e) protein-protein interface (for multimers), or (f) disordered (pLDDT <50). Output: annotated TSV.
11. Snakemake rules: `predict_structure.smk`, `structural_align.smk`.

**Test gate 3:** Predict neprosin (*N. × ventrata*) from sequence using Chai-1. Compute TM-score against 7ZVA (1.80 Å crystal). **Must achieve TM-score >0.85** for the catalytic domain (residues ~130-380 of mature form). If <0.85, predictions for other families without experimental structures cannot be trusted. Also compute per-residue Cα RMSD; mean RMSD should be <2.0 Å for ordered regions.

### Phase 4: FoldX Thermodynamic Analysis

**Scripts to build:**
12. `workflow/scripts/run_foldx_repair.py` — Runs `FoldX --command=RepairPDB` on each predicted structure (5 rounds, as recommended by FoldX documentation). Outputs repaired PDB.
13. `workflow/scripts/run_foldx_scan.py` — For each convergent position identified in Phase 2: mutate from the ancestral amino acid to the derived (carnivorous) amino acid using `FoldX --command=PositionScan`. Also scan the reverse direction (carnivorous → ancestral) to detect asymmetric effects. Temperature set to 298K (standard). Output: raw FoldX output.
14. `workflow/scripts/parse_foldx.py` — Parses FoldX PositionScan output. Extracts total ΔΔG per mutation. Classifies: destabilizing (ΔΔG > 1.0 kcal/mol), neutral (−0.5 to 1.0), stabilizing (ΔΔG < −0.5). These thresholds follow Schymkowitz et al. 2005 recommendations. Output: ΔΔG matrix as TSV.

**Test gate 4:** Run FoldX on the neprosin 7ZVA crystal structure AND the Chai-1 prediction. Compute Pearson correlation of ΔΔG values for all 20 amino acid substitutions at 10 randomly selected positions. Correlation between crystal-derived and prediction-derived ΔΔG should be r > 0.6. If r < 0.6, the predictions are too inaccurate for FoldX analysis and the structure prediction step needs troubleshooting.

### Phase 5: Molecular Docking

**Critical constraint:** AutoDock Vina handles rigid receptors with semiflexible ligands. Peptide substrates must be ≤15 residues. The gliadin 33-mer is too large; use the immunodominant octapeptide fragment PQPQLPYP instead.

**Scripts to build:**
15. `workflow/scripts/prepare_docking.py` — (a) Convert predicted PDB to PDBQT using ADFR Suite `prepare_receptor`. (b) Define docking box centered on active site (identified in Phase 3). Box size = active site bounding box + 10 Å padding on each side. (c) Convert substrate PDB/SDF to PDBQT using Meeko.
16. `workflow/scripts/run_docking.py` — Runs `vina --receptor X.pdbqt --ligand Y.pdbqt --center_x ... --size_x ... --exhaustiveness 32 --num_modes 20`. Per receptor × substrate pair.
17. `workflow/scripts/parse_docking.py` — Extracts top binding affinity (kcal/mol) and top-1 pose RMSD (if reference pose exists for neprosin from 7ZVA).

**Substrates:**

| Substrate | For Enzyme Family | Size | Source |
|-----------|------------------|------|--------|
| PQPQLPYP (gliadin octapeptide) | Neprosins | 8 residues | Immunodominant epitope from gliadin 33-mer |
| GlcNAc₄ (chitin tetramer) | Chitinases | 4 sugar units | Standard chitinase substrate |
| ApA (RNA dinucleotide) | RNase T2 | 2 nt | Minimal substrate for RNase T2 |
| pNPP (p-nitrophenyl phosphate) | Purple acid phosphatases | Small molecule | Standard PAP substrate |
| Hemoglobin octapeptide (VHLTPEEK) | Nepenthesins | 8 residues | Common pepsin-like substrate fragment |

**Test gate 5:** Dock PQPQLPYP into neprosin 7ZVC active site. The predicted pose should place the proline-containing portion within the active-site cleft between the catalytic glutamates (E188, E297 in 7ZVC numbering). Visual inspection in PyMOL/ChimeraX required.

### Phase 6: Electrostatics

**Scripts to build:**
18. `workflow/scripts/run_electrostatics.py` — For each predicted structure: (a) PDB2PQR with PROPKA at pH 2.5, 3.5, and 5.0 (representing the pH range of pitcher fluids across species). (b) APBS to compute electrostatic potential. Output: .dx grid files.

**Analysis:** Compare electrostatic surfaces across orthologs within each family. Question: Do enzymes from species with more acidic pitcher fluid (e.g., *N. rafflesiana*, pH ~2) show systematically different charge distributions than enzymes from less acidic species (e.g., *N. ampullaria*, pH ~5)?

### Phase 7: Expression (LIMITED)

**Only 3 datasets, only target genes.**

| Species | SRA Project | Reference Genome | Citation |
|---------|------------|-----------------|----------|
| *N. gracilis* | PRJDB8591 | GCA_030504385.1 | Saul et al. 2023, *Nature Plants* |
| *Cephalotus follicularis* | PRJDB4470 | GCA_001941015.1 | Fukushima et al. 2017 |
| *Dionaea muscipula* | PRJEB12493 | GCA_032404045.1 | Bemm et al. 2016, *Genome Research* |

**Scripts to build:**
19. `workflow/scripts/quantify_expression.py` — Build Salmon index per species genome. Map reads. Extract TPM for target enzyme genes only. Output: gene × tissue × condition TPM matrix.

### Phase 8: Integration

**Scripts to build:**
20. `workflow/scripts/build_atlas.py` — Creates SQLite database with tables: `enzymes` (metadata), `structures` (PDB paths + pLDDT), `convergence` (positions + annotations), `foldx` (ΔΔG values), `docking` (affinities + poses), `electrostatics` (summary stats), `expression` (TPM matrix).
21. `workflow/scripts/generate_figures.py` — Produces publication figures using matplotlib:
   - Fig 1: Per-family phylogenies with convergent sites highlighted
   - Fig 2: Structural superposition colored by convergent positions (one panel per family)
   - Fig 3: FoldX ΔΔG heatmap (positions × lineages)
   - Fig 4: Docking affinity comparison across orthologs
   - Fig 5: Electrostatic surfaces at pH 2.5 vs 5.0
   - Fig 6: Expression heatmap (if data quality permits)

---

## 8. Script Template

Every script in `workflow/scripts/` must follow this pattern. No exceptions.

```python
#!/usr/bin/env python3
"""One-line description of what this script does.

Reads X from Y, computes Z, writes output to W.
"""

import logging
import sys
from pathlib import Path

import click
import pandas as pd

logger = logging.getLogger(__name__)


@click.command()
@click.option("--input", "-i", type=click.Path(exists=True), required=True,
              help="Path to input file (FASTA/TSV/PDB).")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Path to output file.")
@click.option("--param", type=float, default=1.0,
              help="Description of parameter with units.")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def main(input: str, output: str, param: float, verbose: bool) -> None:
    """Docstring for the CLI command."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )
    input_path = Path(input)
    output_path = Path(output)

    logger.info("Reading input from %s", input_path)
    # ... actual work ...

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # ... write output ...
    logger.info("Wrote output to %s", output_path)


if __name__ == "__main__":
    main()
```

**Rules:**
- All logging goes to stderr.
- All data output goes to files (never stdout in pipelines).
- Every CLI option has a `help` string.
- Paths use `pathlib.Path`, not string concatenation.
- No hardcoded paths. Everything comes from CLI args or config.yaml.
- Type hints on function signatures.
- One function per task; `main()` is thin orchestration.

---

## 9. Snakemake Rule Template

```python
rule example_step:
    """Brief description of this rule."""
    input:
        fasta="results/alignments/{family}.trimmed.fa",
        tree="results/phylogenies/{family}.treefile",
    output:
        tsv="results/convergence/{family}.convergent_sites.tsv",
    params:
        threshold=config["convergence"]["posterior_threshold"],
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/convergence/{family}.log",
    threads: 1
    shell:
        """
        python workflow/scripts/detect_convergence.py \
            --alignment {input.fasta} \
            --tree {input.tree} \
            --output {output.tsv} \
            --threshold {params.threshold} \
            --verbose \
            2> {log}
        """
```

**Rules:**
- Every rule has `log:` capturing stderr.
- Every rule that uses Python packages has `conda:` pointing to the right env.
- Wildcards are `{family}` (enzyme family ID) or `{species}` (species shortname).
- No shell commands longer than 5 lines; complex logic goes in scripts.

---

## 10. Sources

All references cited in this document, verified as of April 2026.

### Core Carnivorous Plant Biology

- Fukushima K, Fang X, Alvarez-Ponce D, et al. (2017) Genome of the pitcher plant *Cephalotus* reveals genetic changes associated with carnivory. *Nature Ecology & Evolution* 1:0059. doi:10.1038/s41559-016-0059
- Wakatake T, Fukushima K. (2025) Transcriptomic prey-capture responses in convergently evolved carnivorous pitcher plants. *New Phytologist* 249:2559–2573. doi:10.1111/nph.70848
- Fukushima K, Pollock DD. (2023) Detecting macroevolutionary genotype–phenotype associations using error-corrected rates of protein convergence. *Nature Ecology & Evolution*. doi:10.1038/s41559-022-01932-7
- Albert VA, Avila-Robledillo L, Fleck S, et al. (2025) Complexity and Innovation in Carnivorous Plant Genomes. NSF PAR 10636585. [Review identifying the outstanding question this project addresses]
- Saul F, Scharmann M, Wakatake T, et al. (2023) Subgenome dominance shapes novel gene space in the Nepenthes pitcher plant. *Nature Plants*. [*N. gracilis* genome]
- Freund M, Graus D, Fleischmann A, et al. (2022) The digestive systems of carnivorous plants. *Plant Physiology* 190:44–59. doi:10.1093/plphys/kiac232
- Pavlovič A. (2025) How the diversity in digestion in carnivorous plants may have evolved. *New Phytologist*. doi:10.1111/nph.70229
- Montero H, Freund M, Fukushima K. (2025) Convergent losses of arbuscular mycorrhizal symbiosis in carnivorous plants. *New Phytologist* 248:2040–2051.

### Enzyme Structures and Biochemistry

- Del Amo-Maestro L, Mendes SR, Rodriguez-Banqueri A, et al. (2022) Molecular and in vivo studies of a glutamate-class prolyl-endopeptidase for coeliac disease therapy. *Nature Communications* 13:4446. doi:10.1038/s41467-022-32215-1 [PDB: 7ZVA, 7ZU8, 7ZVB, 7ZVC]
- Kadek A, Kastanek P, Bezouska K, et al. (2016) Crystallization of nepenthesin I using a low-pH crystallization screen. *Acta Crystallographica Section F* 72:24–28. [Crystallized at 2.8 Å; coordinates NOT deposited in PDB]
- Athauda SBP, Matsumoto K, Rajapakshe S, et al. (2004) Enzymic and structural characterization of nepenthesin, a unique member of a novel subfamily of aspartic proteinases. *Biochemical Journal* 381:295–306.
- Butts CT, Bierma JC, Martin RW. (2016) Novel proteases from the genome of *Drosera capensis*: Structural prediction and comparative analysis. *Proteins* 84:1517–1533.
- Ting TY, Lee WJ, Ramzi AB, Goh HH. (2026) Bioengineered *Saccharomyces cerevisiae* with neprosin for gluten detoxification. *Journal of Future Foods* 6:1227–1238.
- Wall C, Hause F, Grimm W, et al. (2026) A Native Nepenthesin Reactor for Improved Proteolytic Digestion of Intrinsically Disordered Proteins. *ChemBioChem* 27(6).
- Schräder CU, Lee L, Rey M, et al. (2017) Neprosin, a Selective Prolyl Endoprotease for Bottom-up Proteomics and Histone Mapping. *Molecular & Cellular Proteomics* 16:1162–1171.
- Tiew YK, Goh HH. (2022) Neprosin belongs to a new family of glutamic peptidase based on in silico evidence. *Plant Physiology and Biochemistry* 183:23–38.

### Computational Methods

- Emms DM, Kelly S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. *Genome Biology* 20:238. doi:10.1186/s13059-019-1832-y
- Emms DM, Kelly S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons. *Genome Biology* 16:157.
- Katoh K, Standley DM. (2013) MAFFT multiple sequence alignment software version 7. *Molecular Biology and Evolution* 30:772–780.
- Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. (2009) trimAl: a tool for automated alignment trimming. *Bioinformatics* 25:1972–1973.
- Minh BQ, Schmidt HA, Chernomor O, et al. (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference. *Molecular Biology and Evolution* 37:1530–1534.
- Kalyaanamoorthy S, Minh BQ, Wong TKF, et al. (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14:587–589.
- Chai Discovery. (2024) Chai-1: Decoding the molecular interactions of life. bioRxiv 2024.10.10.615955. [Structure prediction]
- Abramson J, Adler J, Dunger J, et al. (2024) Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493–500.
- Elliott LG, Simpkin AJ, Rigden DJ. (2025) ABCFold: easier running and comparison of AlphaFold 3, Boltz-1, and Chai-1. *Bioinformatics Advances* 5:vbaf153.
- van Kempen M, Kim SS, Tumescheit C, et al. (2024) Fast and accurate protein structure search with Foldseek. *Nature Biotechnology* 42:243–246.
- Zhang Y, Skolnick J. (2005) TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Research* 33:2302–2309.
- Schymkowitz J, Borg J, Stricher F, et al. (2005) The FoldX web server: an online force field. *Nucleic Acids Research* 33:W382–W388.
- Guerois R, Nielsen JE, Serrano L. (2002) Predicting changes in the stability of proteins and protein complexes: a study of more than 1000 mutations. *Journal of Molecular Biology* 320:369–387.
- Buß O, Rudat J, Einloft K. (2018) FoldX as Protein Engineering Tool: Better Than Random Based Approaches? *Computational and Structural Biotechnology Journal* 16:25–33. [FoldX benchmarks: 29% success for stabilizing, 69% for destabilizing]
- Eberhardt J, Santos-Martins D, Tillack AF, Forli S. (2021) AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *Journal of Chemical Information and Modeling* 61:3891–3898.
- Trott O, Olson AJ. (2010) AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function. *Journal of Computational Chemistry* 31:455–461.
- Baker NA, Sept D, Joseph S, et al. (2001) Electrostatics of nanosystems. *Proceedings of the National Academy of Sciences* 98:10037–10041.
- Dolinsky TJ, Nielsen JE, McCammon JA, Baker NA. (2004) PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann electrostatics calculations. *Nucleic Acids Research* 32:W665–W667.
- Patro R, Duggal G, Love MI, et al. (2017) Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods* 14:417–419.
- Xue LC, Rodrigues JP, Kastritis PL, et al. (2016) PRODIGY: a web server for predicting the binding affinity of protein-protein complexes. *Bioinformatics* 32:3676–3678.
- Mölder F, Jablonski KP, Letcher B, et al. (2021) Sustainable data analysis with Snakemake. *F1000Research* 10:33.

### Methodological Precedent (FoldX + convergent substitution approach on a different system)

- [Authors]. (2025) Convergence, stability, and thermal adaptation in the rubisco enzyme in plants. bioRxiv 2025.10.08.681247. [Used identical FoldX + convergent substitution pipeline on Rubisco across four plant clades]

---

## 11. What NOT to Build

- **No custom BLAST parser.** OrthoFinder handles all BLAST/DIAMOND internally.
- **No GUI.** Snakemake + Streamlit covers orchestration and visualization.
- **No database ORM.** Raw SQLite with `sqlite3` module. The dataset is <10,000 rows total.
- **No custom file format.** All outputs are TSV, PDB, or JSON.
- **No wrapper scripts around tools that already have Python APIs.** Vina has Python bindings (`from vina import Vina`). FoldX accepts command files directly. Do not re-implement their CLIs.
- **No molecular dynamics.** MD is out of scope. FoldX ΔΔG is sufficient for classifying convergent substitution effects at atlas scale.
- **No AlphaFold2.** Chai-1 or AF3 only. AF2 does not support multimer prediction with the same architecture and is superseded.
- **No homology modeling.** With Chai-1/AF3 available, there is no reason to use SWISS-MODEL or MODELLER. Ab initio prediction from sequence is strictly superior for targets without close templates.

---

## 12. Known Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Chai-1 predictions have low confidence for disordered regions | FoldX results meaningless for those residues | Exclude positions with pLDDT <50 from all downstream analysis |
| FoldX accuracy is moderate (~69% for destabilizing) | Some false positives/negatives in classification | Report distributions, not individual predictions. Focus on aggregate trends across families. Use crystal structure (7ZVA) as internal benchmark. |
| Vina docking scores are not absolute binding affinities | Cannot claim quantitative Kd values | Report only relative rankings within each enzyme family. Compare orthologs to each other, not across families. |
| AF3 weights may not be accessible on Hellbender | Cannot use AF3 | Chai-1 is the primary predictor (pip-installable, Apache 2.0). AF3 is optional. |
| Some carnivorous plant transcriptomes are de novo assemblies (no genome) | Ortholog assignment may include fragmented contigs | Filter by minimum sequence length (>80% of expected full-length) before orthology. |
| Ancestral reconstruction accuracy decreases with phylogenetic distance | Convergent site detection may have false positives at deep nodes | Use posterior probability threshold ≥0.8 (stricter than Fukushima 2017's 0.5) and validate against published sites. |
# ADDENDUM: Methods Review and Corrections

> **Reviewer note to self and to Claude Code:**  
> This project must not become a confirmation loop. Every method choice must be
> challenged against alternatives. If the PI suggests a tool, search for evidence
> that it is suboptimal before agreeing. Do not reinforce ideas — evaluate them.
> This addendum documents a critical re-evaluation of the original CLAUDE.md.

---

## A. Critical Finding: FoldX Alone Is Insufficient

### The Problem

The original guide used FoldX as the sole method for quantifying the effects of convergent substitutions. This is a methodological gap. FoldX answers ONE question: "Does this mutation change the structural stability (ΔΔG) of the folded protein?" It does not answer:

- Is the substitution favored or disfavored by the evolutionary constraints on this protein family?
- Does the substitution alter coevolutionary coupling patterns with other residues?
- Is the substitution's effect on fitness mediated by stability, function, or epistasis?

FoldX has a documented 29% success rate for predicting stabilizing mutations and 69% for destabilizing ones (Buß et al. 2018). In a direct head-to-head, EVmutation outperformed FoldX at identifying variants with enhanced enzyme activity and thermostability in an endoglucanase engineering study (experimental validation of top 20 variants from each method). EVmutation identified variants with up to 3.6-fold activity increase and ΔTm of +2.8°C; FoldX did not.

### What DCA / EVmutation Actually Does (And Why It Is Not a Replacement)

DCA (Direct Coupling Analysis) and its downstream application EVmutation (Hopf et al. 2017, *Nature Biotechnology* 35:128) fit a Potts model to a multiple sequence alignment and infer direct pairwise couplings between residue positions. EVmutation then uses the inferred Hamiltonian to compute a ΔΔE (change in statistical energy) for any single or multiple mutation relative to a reference sequence. This ΔΔE has been shown to correlate with experimentally measured fitness effects across 34 deep mutational scanning datasets (Hopf et al. 2017).

**DCA/EVmutation measures evolutionary constraint. FoldX measures structural stability. These are different quantities.**

A convergent substitution that is:
- FoldX-neutral, EVmutation-favorable → Benefit is functional (specificity, catalysis), not stability
- FoldX-destabilizing, EVmutation-favorable → Stability-function tradeoff under positive selection
- FoldX-neutral, EVmutation-neutral → Likely neutral drift, not adaptive convergence
- FoldX-destabilizing, EVmutation-unfavorable → Likely misidentified as convergent (false positive)

**The comparison between FoldX ΔΔG and EVmutation ΔΔE on the same set of convergent substitutions is itself a novel analysis.** No published study on carnivorous plant enzymes has done this. The Rubisco FoldX study (bioRxiv 2025) used FoldX alone. Butts et al. 2016 (*Drosera* proteases) used docking alone. This dual approach would differentiate the CarnivorEnzyme project from both.

### The MSA Depth Problem

EVmutation requires deep MSAs (typically >1000 effective sequences) for reliable coupling inference. For well-represented families like chitinases (GH19) and aspartic proteases, this is achievable because these are ancient, ubiquitous plant gene families. For neprosins (DUF239), the family is narrower and MSA depth must be verified before running EVcouplings. If Neff < 500, EVmutation predictions will be unreliable for that family.

**Action required per family:** Compute Neff (number of effective sequences at 80% identity threshold) from the alignment. If Neff ≥ 500, run EVcouplings. If Neff < 500, report FoldX only and note the limitation.

---

## B. Corrected Method Stack

### Updated Mutation Effect Prediction (replaces FoldX-only approach)

| Layer | Tool | What It Measures | Requirements | Output |
|-------|------|-----------------|-------------|--------|
| **1. Structural stability** | FoldX 5.0 PositionScan | ΔΔG (kcal/mol) from empirical force field on 3D structure | Predicted or experimental PDB | Per-mutation ΔΔG classification |
| **2. Evolutionary constraint** | EVcouplings/EVmutation (plmc + Python package) | ΔΔE (statistical energy change) from Potts model on MSA | Deep MSA (Neff ≥ 500) | Per-mutation ΔΔE score |
| **3. Zero-shot fitness** | ESM-1v (optional, if compute permits) | Log-likelihood ratio from protein language model | Sequence only (no MSA, no structure) | Per-mutation fitness score |
| **4. Cross-comparison** | Custom analysis script | Scatter plot FoldX ΔΔG vs EVmutation ΔΔE for convergent positions; classify into quadrants | Outputs from layers 1 + 2 | Quadrant classification per convergent site |

**This is the core analytical contribution.** The question "what do convergent substitutions do?" has four possible answers depending on the quadrant:

```
                    EVmutation ΔΔE
                  favorable (+) | unfavorable (-)
                 ________________|________________
FoldX      |                     |                  |
ΔΔG        |   FUNCTIONAL GAIN   |  FALSE POSITIVE  |
neutral    |   (stability ok,    |  (neither stable |
(0)        |    fitness up)       |   nor fit)       |
           |_____________________|__________________|
FoldX      |                     |                  |
ΔΔG        | TRADEOFF / ADAPTIVE | DELETERIOUS      |
destab.    | (stability cost,    | (both structural |
(+)        |  fitness benefit)   |  and evolutionary |
           |                     |  cost)            |
           |_____________________|__________________|
```

### New Scripts Required

Add to `workflow/scripts/`:

- `run_evcouplings.py` — Wrapper for EVcouplings pipeline. Input: per-family alignment. Output: plmc model file (.params) + per-position coupling scores.
- `run_evmutation.py` — Compute ΔΔE for each convergent substitution using the inferred Potts model. Input: .params file + convergent positions TSV. Output: ΔΔE scores per mutation.
- `compare_foldx_evmutation.py` — Merge FoldX ΔΔG and EVmutation ΔΔE. Compute Pearson/Spearman correlation. Classify into quadrants. Generate scatter plots and summary statistics.

### New Snakemake Rules

Add `workflow/rules/evmutation.smk`:

```python
rule run_evcouplings:
    """Infer Potts model from per-family alignment using plmc."""
    input:
        alignment="results/alignments/{family}.trimmed.fa",
    output:
        params="results/evcouplings/{family}.model_params",
        scores="results/evcouplings/{family}.coupling_scores.tsv",
    log:
        "logs/evcouplings/{family}.log",
    conda:
        "../envs/evcouplings.yaml"
    threads: 4
    shell:
        """
        python workflow/scripts/run_evcouplings.py \
            --alignment {input.alignment} \
            --output-params {output.params} \
            --output-scores {output.scores} \
            --threads {threads} \
            --verbose \
            2> {log}
        """

rule run_evmutation:
    """Score convergent substitutions using EVmutation ΔΔE."""
    input:
        params="results/evcouplings/{family}.model_params",
        convergent="results/convergence/{family}.convergent_sites.tsv",
    output:
        tsv="results/evcouplings/{family}.evmutation_scores.tsv",
    log:
        "logs/evmutation/{family}.log",
    conda:
        "../envs/evcouplings.yaml"
    shell:
        """
        python workflow/scripts/run_evmutation.py \
            --model {input.params} \
            --mutations {input.convergent} \
            --output {output.tsv} \
            --verbose \
            2> {log}
        """

rule compare_foldx_evmutation:
    """Cross-compare structure-based and evolution-based mutation effects."""
    input:
        foldx="results/foldx/{family}.ddg_matrix.tsv",
        evmut="results/evcouplings/{family}.evmutation_scores.tsv",
    output:
        tsv="results/integration/{family}.foldx_vs_evmut.tsv",
        plot="results/atlas/figures/{family}_quadrant.pdf",
    log:
        "logs/integration/{family}_compare.log",
    shell:
        """
        python workflow/scripts/compare_foldx_evmutation.py \
            --foldx {input.foldx} \
            --evmutation {input.evmut} \
            --output-tsv {output.tsv} \
            --output-plot {output.plot} \
            --verbose \
            2> {log}
        """
```

### New Conda Environment

Add `workflow/envs/evcouplings.yaml`:

```yaml
name: evcouplings
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - plmc                          # C implementation of pseudolikelihood DCA
  - hmmer=3.4                     # jackhmmer for MSA generation
  - pip
  - pip:
    - evcouplings                 # EVcouplings Python framework
```

---

## C. AF3 Is Available — Corrected Predictor Priority

The user has confirmed local AF3 installation on Hellbender HPC. This changes the recommendation.

**Updated predictor priority:**

1. **AlphaFold3** (primary) — Local installation on Hellbender. Gold standard. Better validated on diverse protein types than Chai-1. Use AF3 JSON input format, 5 seeds per Tier 1 target, 1 per Tier 2.
2. **Chai-1** (fallback) — pip-installable, Apache 2.0. Use if AF3 queue times are prohibitive, or for quick iteration during development.
3. **ABCFold** (comparison) — Run both AF3 and Chai-1 on a subset (all Tier 1 neprosin orthologs) to quantify predictor agreement. Report inter-predictor TM-score.

**Update to config.yaml:**

```yaml
structure:
  primary_predictor: alphafold3     # CHANGED from chai1
  fallback_predictor: chai1
  af3_binary: "/path/to/hellbender/alphafold3"  # Set to actual module path
  n_seeds: 5
  comparison_subset: "tier1_neprosins"  # Run both predictors on this subset
```

**Update to predict_structure.smk:**

The rule must generate AF3-format JSON inputs (not Chai-1 FASTA), run `run_alphafold.py` via SLURM, and extract outputs. The Chai-1 rule should be a separate rule gated by a config flag.

---

## D. Similar Published Projects (Honest Assessment of Novelty)

Projects with overlapping methodology, found via literature search. Evaluate whether CarnivorEnzyme would be scooped or is genuinely novel.

### Directly relevant

1. **Butts, Bierma & Martin 2016** (*Proteins* 84:1517) — Predicted structures of *Drosera capensis* droserasins and nepenthesins using homology modeling (pre-AlphaFold). Performed substrate docking and affinity weight analysis. **Overlap:** Docking of aspartic proteases. **Gap they left:** No cross-lineage comparison, no FoldX, no convergent substitution analysis, no evolutionary constraint analysis.

2. **Fukushima et al. 2017** (*Nat. Ecol. Evol.* 1:0059) — Identified convergent amino acid substitutions in digestive enzymes across carnivorous lineages. Mapped to homology-modeled structures. Found convergent positions are surface-exposed. **Overlap:** Convergent substitution detection + structural mapping. **Gap they left:** No thermodynamic analysis, no docking, no cross-lineage structural comparison using modern predictors, no evolutionary constraint quantification.

3. **Rubisco convergence study** (bioRxiv 2025.10.08.681247) — Used FoldX on spinach RbcL crystal structure mutated to 174 unique sequences from four plant clades. Correlated predicted ΔΔG with growing-season temperature. **Overlap:** FoldX PositionScan on convergent substitutions across plant lineages. **Gap vs CarnivorEnzyme:** Single enzyme (RbcL) vs atlas of 10 families; no DCA/EVmutation; no docking; no electrostatics.

4. **Kim et al. 2025** (*Physiologia Plantarum*) — Used AlphaFold + Foldseek structural phylogeny on Photosystem II D1 subunit diversity. **Overlap:** AlphaFold-based structural evolution of a plant enzyme. **Gap:** Single protein, no mutation effect analysis, no docking.

### Not yet published (potential competition risk)

The Fukushima lab (NIG Japan) is the most likely group to do this work. Their 2025 review explicitly flagged the question this project answers. However, their published work has been genomic/transcriptomic, not structural. They have not published FoldX, docking, or DCA analyses on these enzymes. Risk is moderate — building quickly and posting a preprint is the mitigation.

### Novelty assessment

**No published project combines:**
- Modern structure prediction (AF3/Chai-1) across multiple enzyme families
- Cross-lineage structural comparison using Foldseek
- FoldX thermodynamic analysis of convergent substitutions
- EVmutation evolutionary constraint analysis of the same positions
- FoldX-vs-EVmutation quadrant classification
- Substrate docking comparison across orthologs
- Electrostatic surface comparison correlated with pitcher fluid pH

**The FoldX × EVmutation cross-comparison on convergent sites is the most novel analytical contribution.** No published study has done this on any system, not just carnivorous plants.

---

## E. Updated Validation Strategy

| What | Against | Pass Criteria | Notes |
|------|---------|--------------|-------|
| AF3 neprosin prediction | 7ZVA crystal (1.80 Å) | TM-score > 0.90, mean Cα RMSD < 1.5 Å (ordered regions) | Primary validation |
| Chai-1 neprosin prediction | 7ZVA crystal | TM-score > 0.85 | Secondary; quantify AF3-vs-Chai-1 gap |
| FoldX on AF3 structure vs crystal | ΔΔG Pearson correlation | r > 0.65 | If < 0.65, AF3 structures too inaccurate for FoldX |
| FoldX on crystal vs crystal | ΔΔG self-consistency | r > 0.85 (RepairPDB with different random seeds) | Internal control |
| Convergence detection | Fukushima 2017 known sites (GH19, RNase T2, PAP) | ≥80% recall of published convergent positions | Validates ancestral reconstruction |
| EVcouplings Neff check | Per-family MSA depth | Neff ≥ 500 to proceed with EVmutation | Exclude families below threshold |
| EVmutation sign check | Known deleterious mutations in well-studied homologs | Negative ΔΔE for known pathogenic mutations | Sanity check on model quality |
| Docking gliadin octapeptide into neprosin 7ZVC | Active-site cleft geometry from crystal structure | Proline-containing portion within 4 Å of catalytic glutamates | Visual + distance check |

---

## F. Updated Dependencies

Add to `environment.yml`:

```yaml
  - hmmer=3.4                          # For EVcouplings MSA generation
  - pip:
    - evcouplings                      # DCA + EVmutation pipeline
    - esm                              # ESM-1v zero-shot predictions (optional)
```

Add to external tools:

| Tool | How | Notes |
|------|-----|-------|
| plmc | `conda install -c bioconda plmc` or compile from github.com/debbiemarkslab/plmc | C implementation of pseudolikelihood DCA. Required by EVcouplings. |

---

## G. Updated Figure List

The manuscript should include a FoldX-vs-EVmutation comparison figure. Updated list:

1. Per-family phylogenies with convergent sites highlighted on branches
2. Structural superposition per family, convergent sites colored by surface/buried/active
3. FoldX ΔΔG heatmap (convergent positions × lineages)
4. **EVmutation ΔΔE heatmap (same positions × lineages)**
5. **FoldX ΔΔG vs EVmutation ΔΔE quadrant scatter plots per family** ← NEW, main analytical figure
6. Docking affinity comparison across orthologs
7. Electrostatic surfaces at pH 2.5 vs 5.0
8. Expression heatmap (if data quality permits)

---

## H. Sources for This Addendum

- Hopf TA, Ingraham JB, Poelwijk FJ, et al. (2017) Mutation effects predicted from sequence co-variation. *Nature Biotechnology* 35:128–135. doi:10.1038/nbt.3769
- Morcos F, Pagnani A, Lunt B, et al. (2011) Direct-coupling analysis of residue coevolution captures native contacts across many protein families. *PNAS* 108:E1293–E1301.
- Hopf TA, Green AG, Schubert B, et al. (2019) The EVcouplings Python framework for coevolutionary sequence analysis. *Bioinformatics* 35:1582–1584.
- Buß O, Rudat J, Einloft K. (2018) FoldX as Protein Engineering Tool: Better Than Random Based Approaches? *Comput. Struct. Biotechnol. J.* 16:25–33.
- Livesey BJ, Marsh JA. (2023) Updated benchmarking of variant effect predictors using deep mutational scanning. *Mol. Syst. Biol.* 19:e11474. [EVmutation, EVE, ESM-1v head-to-head on 217 DMS assays]
- Rives A, Meier J, Sercu T, et al. (2021) Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. *PNAS* 118:e2016239118. [ESM-1v]

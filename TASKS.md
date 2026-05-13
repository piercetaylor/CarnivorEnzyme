# CarnivorEnzyme — Tasks

> Date: 2026-05-12
> Companion to [STATUS_2026-05-12.md](STATUS_2026-05-12.md) and [LITERATURE_REVIEW_2026-05-12.md](LITERATURE_REVIEW_2026-05-12.md)

Tasks are organized by tier (urgency × impact), not phase. Each task has an explicit owner, dependency chain, and acceptance criteria. Strike through when complete.

---

## Tier 0 — Pivot Corrections (THIS WEEK)

These three are blocking the science narrative. Resolve before any further compute.

### T0.1 — Swap ESM-2 → ProSST + SaProt
**Owner:** Pierce / Claude Code
**Dependency:** none
**Effort:** ~1 day

- [ ] Create [workflow/envs/plm.yaml](workflow/envs/plm.yaml) with `prosst` and `saprot` pip deps + pytorch ≥ 2.1 + pytorch-cuda
- [ ] Rename [workflow/scripts/score_esm2.py](workflow/scripts/score_esm2.py) → `score_prosst.py`
  - Input: convergent_sites.tsv + alignment FASTA + predicted CIF/PDB structure (NEW — ProSST needs structure tokens)
  - Output TSV columns: `family, aln_position, sequence_id, ancestral_aa, derived_aa, log_p_ancestral, log_p_derived, llr, model="ProSST"`
  - Use Foldseek to extract 3Di tokens from AF3/Chai-1 prediction
- [ ] Create `score_saprot.py` mirroring the same I/O contract with `model="SaProt"`
- [ ] Retain `score_esm2.py` as a sanity baseline but mark it Tier 3 in config
- [ ] Rename rule file `workflow/rules/esm2.smk` → `plm_scoring.smk`. Add three rules: `score_prosst`, `score_saprot`, `score_esm2_baseline`.
- [ ] Update [config/config.yaml](config/config.yaml) `esm2:` block → `plm:` with `primary: prosst`, `secondary: saprot`, `baseline: esm2_t33_650M_UR50D`.
- [ ] Validate on neprosin convergent sites: ProSST LLR should be > 0 for ≥ 4/6 known Fukushima 2017 GH19 convergent positions.

### T0.2 — Add three new genomes + polyploidy annotation
**Owner:** Pierce / Claude Code
**Dependency:** none
**Effort:** ~1 day for species.yaml schema; weeks for actual sequence retrieval (BRAKER2 still required for Tarnita 2023)

- [ ] Add to [config/species.yaml](config/species.yaml):
  - `Sarracenia_purpurea`: carnivory_origin=5, genome=`bioRxiv 2025.12.26.696377` (Albert/Lindqvist chromosome-scale)
  - `Pinguicula_gigantea`: carnivory_origin=6, genome=`bioRxiv 2025.04.05.646448`
  - `Nepenthes_mirabilis`: carnivory_origin=1, genome=`PRJEB86749` (Goh 2025 PacBio HiFi)
- [ ] Extend species.yaml schema with `ploidy: int` and `subgenome_assignments: list[str]` fields
  - *N. gracilis*: ploidy=10, subgenome_assignments=["A","B","C","D","E"]
  - *D. capensis*: ploidy=12
  - *D. muscipula*: ploidy=4
- [ ] Update [config/enzyme_families.yaml](config/enzyme_families.yaml) accession entries — add `subgenome:` tag for polyploid species where >1 homeolog is listed
- [ ] Validate: `fetch_sequences.py` should accept the new entries without error (dry-run mode)

### T0.3 — Enforce FEP scope cap (≤20 sites total)
**Owner:** Claude Code
**Dependency:** none
**Effort:** ~1 hour

- [ ] In [workflow/scripts/run_fep.py](workflow/scripts/run_fep.py), enforce `config.fep.max_targets` strictly. Reject runs exceeding this with explicit error.
- [ ] Add gating logic: site must be in `functional_gain` or `stability_function_tradeoff` quadrant AND ProSST LLR > 0 AND neither aa is Proline.
- [ ] Add tier ordering: Tier 1 families fill first. Within Tier 1, rank by ProSST LLR descending.
- [ ] Verify in [config/config.yaml](config/config.yaml): `fep.max_targets: 20` (currently 10 — bump to 20).

---

## Tier 1 — Critical Path to Preprint

These unlock the primary novel output (Fig. 5: FoldX × ProSST quadrant).

### T1.1 — Run Phase 1 + 2 + 2A on Hellbender for all 5 families
**Owner:** Pierce
**Dependency:** T0.1 complete
**Effort:** ~2 days HPC time

- [ ] On Hellbender: rsync the latest config/ + workflow/ from local
- [ ] `snakemake --use-conda phase1 --cores 4 --executor slurm --profile config/slurm/`
- [ ] `snakemake --use-conda phase2 --executor slurm`
- [ ] `snakemake --use-conda phase2A --executor slurm` (convergence detection)
- [ ] Test gate 2: nepenthesin tree topology — N. gracilis Nep I/II monophyletic, Cephalotus outside Caryophyllales clade
- [ ] Test gate 3: GH19 convergence — ≥ 80% recall of Fukushima 2017 Fig. 3a positions

### T1.2 — Implement structure prediction (`predict_af3.py` + `predict_chai1.py`)
**Owner:** Claude Code
**Dependency:** none (parallel with T1.1)
**Effort:** ~2 days

- [ ] [workflow/scripts/predict_af3.py](workflow/scripts/predict_af3.py): two-stage (MSA CPU → inference GPU); 5 seeds Tier 1
- [ ] [workflow/scripts/predict_chai1.py](workflow/scripts/predict_chai1.py): fallback when AF3 queue too deep
- [ ] [workflow/scripts/assess_structure.py](workflow/scripts/assess_structure.py): pLDDT extraction, TM-align vs reference PDB
- [ ] [workflow/rules/predict_structure.smk](workflow/rules/predict_structure.smk): full Snakemake rule with SLURM submission
- [ ] Test gate 4: neprosin AF3 prediction vs 7ZVA — TM-score > 0.90

### T1.3 — Implement FoldX pipeline (pH-aware physics axis)
**Owner:** Claude Code
**Dependency:** T1.2 (needs predicted PDBs)
**Effort:** ~2 days

- [ ] Tighten `config.structure.plddt_exclude: 50 → 70` in [config/config.yaml](config/config.yaml) (per [FOLDX_REVIEW_2026-05-12.md](FOLDX_REVIEW_2026-05-12.md))
- [ ] [workflow/scripts/run_foldx_repair.py](workflow/scripts/run_foldx_repair.py): 5 rounds of RepairPDB
- [ ] [workflow/scripts/run_foldx_scan.py](workflow/scripts/run_foldx_scan.py): PositionScan at pH 2.5, 3.5, 5.0
  - **Upstream:** call PypKa (Reis 2024, *NAR* 52:W294) for protonation at each pH instead of PROPKA — better validated on AF inputs
- [ ] [workflow/scripts/parse_foldx.py](workflow/scripts/parse_foldx.py): ΔΔG matrix → classified TSV
- [ ] [workflow/rules/foldx.smk](workflow/rules/foldx.smk)
- [ ] Test gate 5: ΔΔG on neprosin 7ZVA crystal vs AF3 prediction — Pearson r > 0.65

### T1.3b — Implement deep-learning ΔΔG ensemble (NEW — see FOLDX_ALTERNATIVES_AF3_2026-05-12.md)
**Owner:** Claude Code
**Dependency:** T1.2 (needs predicted PDBs); parallel with T1.3
**Effort:** ~3 days total
**Why:** FoldX 5.1 PCC ~0.40 on Megascale vs ThermoMPNN PCC ~0.75, SPURS ~0.54 on AF-native Human Domainome. Single-method FoldX is no longer reviewer-defensible without DL backstop.

- [ ] Create [workflow/envs/stability_ensemble.yaml](workflow/envs/stability_ensemble.yaml) with `spurs`, `rasp`, `catopt`, `ddgun3d` pip/conda deps
- [ ] `workflow/scripts/score_spurs.py` — **AF-native ΔΔG**, primary DL axis. Cao et al. 2025 *Nat. Commun.* doi:10.1038/s41467-025-67609-4. Code: [github.com/luo-group/SPURS](https://github.com/luo-group/SPURS)
- [ ] `workflow/scripts/score_rasp.py` — saturation-scale stability map per ortholog. Blaabjerg 2023 *eLife* doi:10.7554/eLife.82593
- [ ] `workflow/scripts/score_catopt.py` — **Δ pHopt per mutation** (biological direction of pH adaptation). Wang 2025 *ACS Synth. Biol.* doi:10.1021/acssynbio.5c00679
- [ ] `workflow/scripts/score_ddgun3d.py` — untrained baseline, robust to AF structure quality. Montanucci 2019
- [ ] `workflow/rules/stability_ensemble.smk` — single rule producing 4 stability TSVs per family
- [ ] Test: at least 3 of 4 axes agree on sign of ΔΔG for known stabilizing/destabilizing positions (validate against 7ZVA crystal mutations from literature)

### T1.4 — Implement multi-axis quadrant comparison (PRIMARY NOVEL OUTPUT — revised scope)
**Owner:** Claude Code
**Dependency:** T0.1 + T1.3 + T1.3b
**Effort:** ~1 day

- [ ] Rename `compare_foldx_evmutation.py` → [compare_stability_evolution.py](workflow/scripts/compare_stability_evolution.py)
- [ ] Merge **four stability axes** (FoldX pH 2.5/3.5/5.0, SPURS, RaSP, CatOpt) with **three evolutionary axes** (ProSST, SaProt, EVE/EVcouplings) on (family, position, ancestral_aa, derived_aa)
- [ ] Compute consensus stability score (mean of normalized z-scores across 4 axes)
- [ ] Compute consensus evolutionary score (mean of normalized LLR across 3 axes)
- [ ] Quadrant classification on consensus axes:
  - Q1 (functional_gain): consensus_stab ≤ 1.0 AND consensus_evo > 0
  - Q2 (stability_function_tradeoff): consensus_stab > 1.0 AND consensus_evo > 0
  - Q3 (deleterious): consensus_stab > 1.0 AND consensus_evo ≤ 0
  - Q4 (neutral_drift / false_positive): consensus_stab ≤ 1.0 AND consensus_evo ≤ 0
- [ ] Generate publication-quality scatter plot, 5 panels (one per Tier 1 family) → Fig. 5 of manuscript
- [ ] Add supplementary plot: per-axis disagreement heatmap (flags sites where methods disagree — these are uncertain)

### T1.4 — Implement FoldX × ProSST quadrant comparison (PRIMARY NOVEL OUTPUT)
**Owner:** Claude Code
**Dependency:** T0.1 + T1.3
**Effort:** ~1 day

- [ ] Rename `compare_foldx_evmutation.py` → [compare_foldx_prosst.py](workflow/scripts/compare_foldx_prosst.py)
- [ ] Merge FoldX ΔΔG TSV with ProSST LLR TSV on (family, position, ancestral_aa, derived_aa)
- [ ] Quadrant classification:
  - Q1 (functional_gain): ΔΔG ≤ 1.0 kcal/mol AND LLR > 0
  - Q2 (stability_function_tradeoff): ΔΔG > 1.0 AND LLR > 0
  - Q3 (deleterious): ΔΔG > 1.0 AND LLR ≤ 0
  - Q4 (false_positive): ΔΔG ≤ 1.0 AND LLR ≤ 0 (neutral drift candidates)
- [ ] Generate publication-quality scatter plot (matplotlib, 300 dpi, serif font)
- [ ] One panel per Tier 1 family (5 panels total → Fig. 5 of manuscript)

---

## Tier 2 — Manual Data Retrieval (Per temporary_TODO_42526.md)

These are not blockers for the quadrant analysis (the 84 verified accessions already cover ≥ 3 origins per family) but are needed for full coverage.

### T2.1 — Run BRAKER2 on 5 Tarnita 2023 genomes (Hellbender)
**Owner:** Pierce (manual)
**Dependency:** Hellbender SLURM access
**Effort:** ~30 min setup + 12–24h job runtime
**Detailed steps:** [temporary_TODO_42526.md](temporary_TODO_42526.md) section 1

- [ ] `bash workflow/scripts/braker2/00_download_genomes.sh`
- [ ] `bash workflow/scripts/braker2/01_download_hints.sh`
- [ ] `bash workflow/scripts/braker2/02_submit_all.sh`
- [ ] Wait for jobs (12–24 h)
- [ ] `bash workflow/scripts/braker2/03_extract_proteins.sh`
- [ ] `bash workflow/scripts/braker2/04_run_hmmer_scan.sh`
- [ ] Update [config/enzyme_families.yaml](config/enzyme_families.yaml) TODO entries → `local:` paths

### T2.2 — ENA download of *N. mirabilis* proteins (7 gene IDs)
**Owner:** Pierce (manual)
**Dependency:** internet + samtools
**Effort:** ~1 h
**Detailed steps:** [temporary_TODO_42526.md](temporary_TODO_42526.md) section 2

- [ ] `curl -L "https://www.ebi.ac.uk/ena/browser/api/fasta/PRJEB86749?download=true" -o resources/accessions/n_mirabilis_all_proteins.fasta`
- [ ] `samtools faidx resources/accessions/n_mirabilis_all_proteins.fasta`
- [ ] Extract the 7 gene IDs (jg4646.t1, jg6978.t1, jg11529.t1, jg39601.t1, jg12124.t1, jg35245.t1, jg38588.t1)
- [ ] Add `local:` support to [workflow/scripts/fetch_sequences.py](workflow/scripts/fetch_sequences.py) (Tier 1, T1.5 below)
- [ ] Update YAML TODOs → `local:` paths

### T2.3 — JGI download of *U. gibba* gene models
**Owner:** Pierce (manual; JGI account needed)
**Dependency:** JGI login
**Effort:** ~half day
**Detailed steps:** [temporary_TODO_42526.md](temporary_TODO_42526.md) section 3

- [ ] Manually download `Utrgib1_proteins.aa.fasta.gz` from genome.jgi.doe.gov (JGI portal Utrgib1)
- [ ] Run hmmbuild + hmmsearch for GH19, Metallophos, RNase T2, Asp protease
- [ ] Extract top hits per family
- [ ] Update YAML TODOs → `local:` paths

### T2.4 — CDD domain verification (8 accessions flagged "hypothetical protein")
**Owner:** Pierce (manual)
**Dependency:** none
**Effort:** ~1 h (web tool batch submit)
**Detailed steps:** [temporary_TODO_42526.md](temporary_TODO_42526.md) section 6

- [ ] Submit each accession to https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi
- [ ] Confirm domain assignments and update YAML notes

### T1.5 — Add `local:` prefix support to fetch_sequences.py
**Owner:** Claude Code
**Dependency:** none (but unblocks T2.1, T2.2, T2.3)
**Effort:** ~1 hour

- [ ] Modify [workflow/scripts/fetch_sequences.py](workflow/scripts/fetch_sequences.py) to recognize `local:path/to/file.fasta` accession syntax
- [ ] Behavior: copy or symlink the local file into `results/sequences/{family}/{species}.fa` instead of querying Entrez
- [ ] Validate length filter still applies to local sequences
- [ ] Add test in [tests/test_fetch.py](tests/test_fetch.py)

---

## Tier 3 — Downstream Pipeline Stubs (After Quadrant Lands)

Implement after the Fig. 5 quadrant figure is generated for at least one family. These are not preprint blockers but are needed for the full atlas.

### T3.1 — EVcouplings / EVE for deep-MSA families
**Effort:** ~2 days
- [ ] [workflow/scripts/run_evcouplings.py](workflow/scripts/run_evcouplings.py) with `plmc` + Neff ≥ 500 gate
- [ ] Consider EVE as upgrade per Frazer 2021 — stronger calibration than EVmutation
- [ ] [workflow/scripts/run_evmutation.py](workflow/scripts/run_evmutation.py) → score convergent sites with Potts model
- [ ] Gate output: only families with Neff ≥ 500. For neprosins (low Neff), skip and rely on ProSST/SaProt.

### T3.2 — ACEP embedding-distance convergence (Stiffler 2025 PNAS)
**Effort:** ~1 week
- [ ] New script: `workflow/scripts/run_acep.py`
- [ ] Compute PLM embeddings (ESM-2 or ProSST) for each ortholog
- [ ] Compute pairwise embedding distance across independent carnivorous lineages
- [ ] Test whether convergent lineages are closer in embedding space than expected under null
- [ ] Complement to site-level ωC (Fukushima 2023)

### T3.3 — Substrate docking pipeline
**Effort:** ~3 days
- [ ] [workflow/scripts/prepare_docking.py](workflow/scripts/prepare_docking.py) (ADFR + Meeko)
- [ ] [workflow/scripts/run_docking.py](workflow/scripts/run_docking.py) (Vina Python API)
- [ ] [workflow/scripts/parse_docking.py](workflow/scripts/parse_docking.py)
- [ ] [workflow/rules/docking.smk](workflow/rules/docking.smk)
- [ ] Test gate 6: PQPQLPYP into neprosin 7ZVC — proline portion within 4 Å of E188/E297

### T3.4 — Electrostatics (PDB2PQR + APBS)
**Effort:** ~2 days
- [ ] [workflow/scripts/run_electrostatics.py](workflow/scripts/run_electrostatics.py): pH 2.5 / 3.5 / 5.0
- [ ] [workflow/rules/electrostatics.smk](workflow/rules/electrostatics.smk)
- [ ] Correlation analysis: surface potential at pH 2.5 vs species pitcher_pH

### T3.5 — Expression quantification (Salmon)
**Effort:** ~2 days
- [ ] [workflow/scripts/quantify_expression.py](workflow/scripts/quantify_expression.py)
- [ ] [workflow/rules/expression.smk](workflow/rules/expression.smk)
- [ ] Three datasets only: PRJDB8591 (*N. gracilis*), PRJDB4470 (*Cephalotus*), PRJEB12493 (*Dionaea*)

### T3.6 — Atlas SQLite + 10 figures
**Effort:** ~3 days
- [ ] [workflow/scripts/build_atlas.py](workflow/scripts/build_atlas.py): schema in [AGENTS.md](AGENTS.md) Phase 9
- [ ] [workflow/scripts/generate_figures.py](workflow/scripts/generate_figures.py): all 10 publication figures
- [ ] [workflow/rules/integrate.smk](workflow/rules/integrate.smk)

---

## Tier 4 — Strategic / External

### T4.1 — Wet-lab validation collaboration
**Owner:** Pierce
**Effort:** months from initiation; ~1 figure output
**Decision deadline:** before bioRxiv submission

- [ ] Draft outreach email to Goh lab (PRJEB86749) or Brodelius lab
- [ ] Proposal: express 2–3 ProSST/FoldX-quadrant-winner neprosin variants; measure activity vs gliadin peptide; compare to wild-type
- [ ] Even one experimental validation panel converts paper from 30-cite to 100-cite class

### T4.2 — Manuscript outline
**Owner:** Pierce
**Effort:** ~1 week
**Target:** initiate by 2026-06-01

- [ ] Draft section headings
- [ ] Cite DETANGO (Ding 2026), Chugunov 2025, Rubisco paper (bioRxiv 2025.10.08.681247), Stiffler PNAS 2025 — see [LITERATURE_REVIEW_2026-05-12.md](LITERATURE_REVIEW_2026-05-12.md)
- [ ] Distinguish CarnivorEnzyme: multi-family + cross-lineage + convergence focus

### T4.3 — Tests (currently all stubs)
**Effort:** ~1 day per script
- [ ] Replace placeholder tests in [tests/](tests/) with real fixture-based unit tests
- [ ] Prioritize: `test_convergence.py`, `test_foldx_parse.py`, `test_evmutation.py` (→ rename `test_prosst.py`)

---

## Completion ledger

When a task is complete, change `- [ ]` to `- [x]` and add date.
When a tier is fully done, note here:

- [ ] Tier 0 complete (pivot corrections)
- [ ] Tier 1 complete (preprint critical path)
- [ ] Tier 2 complete (manual data retrieval)
- [ ] Tier 3 complete (full downstream pipeline)
- [ ] Tier 4 complete (manuscript + validation)

# CarnivorEnzyme — TODO

> Date: 2026-08-26
> Supersedes: `TASKS.md` (archived, see `archive/TASKS.md`)
> Companion documents: [PROJECT_STATUS.md](PROJECT_STATUS.md), [PROJECT_PLAN.md](PROJECT_PLAN.md)

Tasks are organized by tier (blocking → strategic → downstream), not by phase number. Every task
carries the same gate format so "done" is checkable, not a judgment call:

- **Owner** — who does this
- **Dependency** — what must be true/complete before starting
- **Effort** — rough estimate
- **Entry condition** — the state that makes this task ready to start
- **Steps** — what to actually do
- **Verification** — the exact command/check to run
- **Pass/fail criteria** — the specific bar that decides done vs. not-done

Check off `- [ ]` → `- [x]` with a date when a step completes. Don't mark a tier's ledger entry
done until every task in it is checked.

---

## Tier 0 — Blocking HPC verification (must clear before any other tier)

**Why this is first:** the ancestral-reconstruction architecture fix (commit `9b1a8c5`,
2026-08-24) that convergence detection now depends on has only been verified against synthetic
data and a locally-downloaded Windows IQ-TREE binary — never against Hellbender or this project's
real alignments. Nothing downstream of Phase 2 is trustworthy until this closes.

### T0.1 — Run the real phylogeny→convergence chain on Hellbender
**Owner:** Pierce
**Dependency:** none
**Effort:** ~1 day HPC time
**Entry condition:** Hellbender SLURM access confirmed working (see Hellbender SSH policy in
project memory — batch commands, don't SSH unprompted)

- [ ] rsync current `config/` + `workflow/` to Hellbender
- [ ] Run `phylogeny.smk` → `convergence.smk` for `chitinases_gh19_class_iv` (best-populated Tier-1
      family)
- [ ] Confirm the IQ-TREE binary name/flags actually match the Hellbender installation (`iqtree` vs
      `iqtree2` vs `iqtree3` — `config.yaml phylogeny.binary` currently assumes `iqtree3`, verified
      only against a local Windows 3.1.3 binary)

**Verification:** `results/convergence/chitinases_gh19_class_iv.convergent_sites.tsv` exists and
is non-empty.
**Pass/fail:** Rows in the output are real (traceable to actual alignment positions and lineages),
not fabricated — cross-check a handful of rows by hand against the alignment. Fabricated/empty
output with no error = fail, re-open `audit/11_ancestral_reconstruction_architecture_fix.md`.

### T0.2 — Fix the `ChildIOException` Snakemake DAG bug
**Owner:** Claude Code
**Dependency:** none (parallel with T0.1)
**Effort:** ~2 hours
**Entry condition:** none

- [x] In `workflow/rules/retrieve.smk`, `fetch_all_sequences`'s output (`results/sequences/`) is an
      ancestor directory of `combine_family_sequences`'s output — restructure one of the two rules'
      output paths so neither is a parent of the other (2026-09-02: `combined.fa` moved to
      `results/family_fasta/{family}.combined.fa`; see `audit/16_snakemake_dag_childio_fix.md`)
- [x] Also fixed in the same pass: `esm2.smk` declared an input `all_sequences.fa` that **no rule
      produces**, which would have broken `phase3a` as soon as the ChildIOException stopped masking
      it. Repointed to the same combined FASTA (2026-09-02)

**Verification:** `snakemake -n` (dry run) targeting any rule past Phase 1 no longer raises
`ChildIOException`. **Done 2026-09-02** — verified by execution against every phase target: the
T0.1 target builds an 8-job DAG in the correct order, `phase2` 26 jobs, `phase3a` 34 jobs;
`phase4a`/`phase5d`/`all` now fail only on stub rules, as expected. 13/13 tests still pass.
**Pass/fail:** Dry run reaches the DAG-construction stage without this specific error. (Later
stages may still fail on stub rules — that's expected until Tier 2 closes, not this task's
concern.)

### T0.3 — Audit `extract_ancestor.py` for the same node-naming bug class
**Owner:** Claude Code
**Dependency:** T0.1 (real data makes this easier to test against)
**Effort:** ~half day
**Entry condition:** T0.1 complete

- [ ] Check `extract_ancestor.py` for the same fallback patterns (levelorder-index naming,
      substring `.state` matching) that were removed from `detect_convergence.py` in the Aug-24 fix
- [ ] If present, apply the same fix pattern: validate node-by-node against the `.asr.state` file,
      raise loudly on mismatch, never invent a node name

**Verification:** Run `extract_ancestor.py` against real Tier-0-verified output; deliberately feed
it a mismatched tree/state pair and confirm it raises rather than silently producing a wrong MRCA
sequence.
**Pass/fail:** No silent-fallback path remains reachable.

**Ledger:** - [ ] Tier 0 complete (all three tasks checked)

---

## Tier 0.5 — Repository hygiene (cheap, parallel-safe, do alongside Tier 0)

### T0.5.1 — Remove the stray temp config file
**Owner:** Claude Code | **Effort:** 5 min
- [ ] Delete `config/enzyme_families.yaml.tmp.41128.1776193383614` after confirming it's not
      referenced anywhere (`grep -r` for the filename) and diffing it against the real
      `enzyme_families.yaml` to confirm nothing unique was lost
**Verification:** file gone, `test_config_consistency.py` still passes.

### T0.5.2 — Resolve `hellbenderModules.txt`
**Owner:** Pierce | **Effort:** 15 min
- [ ] Either populate it with the actual `module load` commands needed on Hellbender (iqtree,
      conda, etc.) or delete it if nothing depends on it
**Verification:** File is either non-empty and accurate, or absent — not a silent 0-byte stub.

### T0.5.3 — Archive superseded status/task docs
**Owner:** Claude Code | **Effort:** 15 min
- [x] Move `STATUS_2026-05-12.md`, `TASKS.md`, `temporary_TODO_42526.md`,
      `PROJECT_RESTRUCTURE_2026-05-12.md` into a new top-level `archive/` folder (2026-08-26)
- [x] `PROGRESS.md` intentionally **not** moved — it's gitignored (`.gitignore:111`) and was never
      tracked; forcing it into version control under `archive/` would fight that existing
      convention. Left at repo root as-is. (2026-08-26)
- [x] Prepend each moved file with a one-line header: "Superseded by PROJECT_STATUS.md /
      PROJECT_PLAN.md / TODO.md, 2026-08-25. Preserved for history, not current." (2026-08-26)
**Verification:** `git status` shows the moves; no broken relative links remain inside the archived
files themselves (links to each other within `archive/` are fine to leave as-is).

### T0.5.4 — Fix README.md's known drift
**Owner:** Claude Code | **Effort:** ~1 hour
**Dependency:** T0.5.3 (so README can point at the new doc set)
- [ ] Fix the stale "nine independent lineage origins" reference (should reflect the corrected
      single-Caryophyllales-origin count per CLAUDE.md §4)
- [ ] Update README's doc-map links to point at PROJECT_STATUS.md / PROJECT_PLAN.md / TODO.md
      instead of the archived docs
- [ ] Reconcile README's phase numbering (currently describes a 9-phase DAG matching the old
      STATUS_2026-05-12.md numbering) with CLAUDE.md §7's phase numbering
**Verification:** grep README.md for "nine independent" and for links to archived filenames —
both should return nothing.

### T0.5.5 — Make `rule all` the default Snakemake target (NEW 2026-09-02)
**Owner:** Claude Code | **Effort:** 15 min
Surfaced while verifying T0.2. `snakemake -n` with no target resolves to **`rule phase1`**, not
`rule all` — Snakemake takes the first rule in the file, and `phase1` is at Snakefile:83 while `all`
is at Snakefile:141. The Snakefile's own usage docstring says `snakemake --use-conda --cores 8`
runs the full pipeline; it does not, it downloads sequences and exits 0. Left unfixed deliberately —
changing the workflow's default target is a behavioural decision, not part of a DAG bugfix.
- [ ] Either move `rule all` above the phase targets, or add `default_target: True` to it
- [ ] Reconcile the Snakefile usage docstring with whichever is chosen
**Verification:** `snakemake -n` with no target builds the full-pipeline DAG (or fails on stub
rules), not a 2-job Phase 1 DAG.

**Ledger:** - [ ] Tier 0.5 complete

---

## Tier 1 — Method-stack migration (2026-08-25 decisions, as revised 2026-09-01)

> **T1.1 and T1.2 were rewritten on 2026-09-01.** The stability half of the 2026-08-25 decision was
> reversed by `audit/15_stability_predictor_audit.md`: FoldX is reinstated as primary with RaSP as a
> cross-check, SPURS is demoted to optional, and the pH claim moves to Phase 6 electrostatics. The
> evolutionary-axis decisions (T1.3–T1.5, ProSST/SaProt/EVE) are unaffected.

**Dependency for the whole tier:** Tier 0 complete (don't migrate tooling on top of an unverified
pipeline core).

### T1.1 — Correct CLAUDE.md's Method Justification Table (REVISED 2026-09-01)
**Owner:** Claude Code | **Effort:** ~2 hours
- [x] FoldX row rewritten: correct citation (Delgado et al. 2025, not the fabricated "Botte et al."),
      corrected numbers (R 0.693 → 0.706, RMSE 1.252 — the 0.711/1.238 figures describe the
      median-of-five-runs protocol, not the force-field revision), and honest limitations
- [x] RaSP row added as the structure-sensitivity cross-check
- [x] CpHMD row corrected: real citations (Aho et al. 2022; Buslaev et al. 2022; Jansen et al. 2024),
      and the false "native / GROMACS 2024.2 / conda-installable" claims removed
- [ ] Replace the ESM-2 row with ProSST (primary) + SaProt (secondary) + ESM-2 (Tier-3 baseline)
- [ ] Replace the EVmutation/EVcouplings row with EVE (Frazer et al. 2021)
- [ ] Update CLAUDE.md §7 Phase 4 and Phase 5E test gates to match PROJECT_PLAN.md §4
- [ ] Rewrite addendum §A/§B — they still argue the FoldX×EVmutation quadrant with kcal/mol
      thresholds; the stability axis is now a within-family rank (PROJECT_PLAN.md §2.5)
- [ ] Revise addendum §D's novelty claim to the tightened version in PROJECT_PLAN.md §5

**Verification:** `grep -n "Botte\|EVmutation\|EVcouplings" CLAUDE.md` — remaining hits only in
explicitly-historical context.
**Pass/fail:** No row in the table cites a fabricated author or presents a superseded tool as active.

### T1.1b — Crossref sweep of every remaining citation (NEW 2026-09-01)
**Owner:** Claude Code | **Effort:** ~half day | **Priority: do before any manuscript text exists**
Five fabricated references have been found so far — "Botte et al." (FoldX), "Cai"/"Cao" (SPURS),
"Ressa" (JanusDDG), and two invented GROMACS CpHMD refs sharing one bogus page number. In every case
the DOI was correct and the author name was invented, which is the signature of citations assembled
from search snippets rather than read. Assume the rest of the corpus has the same defect rate.
- [ ] Extract every citation from `CLAUDE.md`, `README.md`, `PROJECT_PLAN.md` and `literature_review/`
- [ ] Resolve each DOI against Crossref; flag any where the deposited first author does not match
- [ ] Fix or delete every mismatch; record the sweep in `audit/`

**Pass/fail:** Zero unresolved citations remain in the four locations above.

### T1.2 — Complete the FoldX phase and add RaSP (REVISED 2026-09-01)
**Owner:** Claude Code | **Dependency:** T1.1 | **Effort:** ~2 days
> Reversal note: the prior version of this task said "Delete (not complete) `run_foldx_repair.py`,
> `run_foldx_scan.py`, `parse_foldx.py`, `foldx.smk`, `workflow/envs/foldx.yaml`". **Do not delete
> them — complete them.** FoldX 5.1 is licensed and installed.
- [ ] Complete `run_foldx_repair.py` (RepairPDB, 5 rounds) and `run_foldx_scan.py` (PositionScan at
      `config.foldx.reference_ph`), and fill in `foldx.smk` / `workflow/envs/foldx.yaml`
- [ ] `parse_foldx.py` must emit **within-family percentile rank** alongside raw kcal/mol, and must
      carry a per-position relative-SASA column so surface and buried sites are never pooled in a
      figure or a summary statistic
- [ ] Add the supplementary pH sweep behind `config.foldx.supplementary_ph_sweep`, clearly labelled
      as a sensitivity check and excluded from main figures
- [ ] `workflow/envs/stability.yaml` + `workflow/scripts/score_rasp.py` + `workflow/rules/rasp.smk`,
      matching `parse_foldx.py`'s output shape so the two are directly rank-comparable
- [ ] `score_spurs.py` remains optional; build it only after FoldX and RaSP are both landing results
- [ ] Port the audit's scan harness (SASA-stratified position selection + pH comparison) into
      `tests/` as a regression fixture, so the surface/buried degeneracy check is repeatable

**Verification:** reinstated Test Gate 4 — FoldX ΔΔG on the 7ZVA crystal vs. the AF3 neprosin
prediction must correlate at r > 0.6 across matched positions; and FoldX vs. RaSP rank correlation
ρ > 0.5 within a family, with disagreement reported rather than resolved.
**Pass/fail:** as stated. Note the audit already ran FoldX 5.1 end-to-end on 7ZVA chain A (RepairPDB
+ PositionScan, 87 mutations, three pH values), so the tool path itself is known to work.

### T1.2b — Promote Phase 6 electrostatics to the primary pH axis (NEW 2026-09-01)
**Owner:** Claude Code | **Dependency:** T2.1 (needs structures) | **Effort:** ~2 days
The project's pH claim can be carried by neither the FoldX pH sweep (effect ~10× below the tool's
error bar) nor CpHMD as configured (unreleased fork; 3 pH points cannot yield a pKa). Surface
electrostatics is what Fukushima 2017's own interpretation and the acidophile-adaptation literature
both point to.
- [ ] Move `run_electrostatics.py` out of Tier 2 ordering into Tier 1 priority
- [ ] PDB2PQR + PROPKA + APBS at pH 2.5 vs 5.0; report surface-charge distribution and
      active-site electrostatic potential per ortholog
- [ ] State PROPKA's real error in Methods: 0.79 pH units for Asp/Glu, **1.37 for buried carboxylates**
- [ ] Do not adopt PypKa — it beats PROPKA by 0.04 pH units and beats doing nothing by 0.02

**Verification:** surface potential differences are computable and directionally interpretable for
≥1 Tier-1 family, and are reported against the *measured* pitcher-pH subset only.

### T1.2c — Fix the phenotype table or downgrade the comparative claim (NEW 2026-09-01)
**Owner:** Pierce | **Effort:** ~half day
`config/species.yaml` has `pitcher_pH` for **9 of 30 species**, spanning 2.5–5.5, in 5
`carnivory_origin` groups with four of the nine inside *Nepenthes*. A PGLS against pitcher pH is
n ≈ 9 with heavy phylogenetic clustering.
- [ ] Either expand the phenotype table from published sources (Adlassnig et al. 2011 is the
      obvious candidate), or
- [ ] Report the comparative analysis as descriptive, with no headline p-value

### T1.2d — Carried-over gaps from the archived paper-structure doc (NEW 2026-09-01)
`archive/ENZYME_EVOLUTION_PAPER_STRUCTURE.md` identified two analyses that reviewers in this field
expect and that the project does not have. They survive the archival.
- [ ] **dN/dS selection tests** — PAML codeml branch-site, or HyPhy aBSREL/BUSTED/RELAX. Requires a
      nucleotide CDS retrieval path that does not currently exist
- [ ] **PGLS phenotype correlation** — gated on T1.2c

### T1.3 — Build `score_prosst.py` + `score_saprot.py`
**Owner:** Claude Code | **Dependency:** T1.1; hard-blocked on Phase 3 (Tier 2) for real execution, but the script itself can be written now | **Effort:** ~1 day
- [ ] `workflow/envs/plm.yaml` — prosst + saprot + pytorch ≥2.1 + pytorch-cuda
- [ ] `score_prosst.py` — needs Foldseek 3Di tokens extracted from a predicted structure (new
      dependency: Foldseek must run before this script, not just before Phase 3's own structural
      alignment step)
- [ ] `score_saprot.py` mirroring the same I/O contract
- [ ] Rename `workflow/rules/esm2.smk` → `plm_scoring.smk`; add `score_prosst`, `score_saprot`,
      `score_esm2_baseline` rules
- [ ] Update `config.yaml`: rename `esm2:` block → `plm:` with `primary: prosst`,
      `secondary: saprot`, `baseline: esm2_t33_650M_UR50D`

**Verification:** Cannot run for real until Tier 2's structure-prediction phase produces an actual
predicted structure. Write and unit-test the script logic now; defer the end-to-end run.
**Pass/fail (deferred to Tier 2):** ProSST LLR positive for ≥4/6 known Fukushima 2017 GH19
convergent positions (PROJECT_PLAN.md §4).

### T1.4 — Rename/rebuild the DCA scripts for EVE
**Owner:** Claude Code | **Dependency:** T1.1 | **Effort:** ~1 day
- [ ] Rename `run_evcouplings.py` → an EVE-targeting equivalent (or rewrite in place — EVE's model-
      fitting step differs from plmc/EVcouplings, check before assuming a pure rename works)
- [ ] Rename `run_evmutation.py` → EVE scoring equivalent
- [ ] Rename `evmutation.smk` → `eve.smk`
- [ ] Update `workflow/envs/evcouplings.yaml` → `eve.yaml` with the correct EVE deps
- [ ] `config.yaml`: rename `evcouplings:` block → `eve:`, keep the `min_neff: 500` gate unchanged

**Verification:** Neff computation still gates correctly (same logic, new tool name).
**Pass/fail:** Negative ΔΔE for a known deleterious mutation in a well-studied homolog
(PROJECT_PLAN.md §4).

### T1.5 — Rebuild the quadrant-comparison script
**Owner:** Claude Code | **Dependency:** T1.2, T1.3, T1.4 | **Effort:** ~1 day
- [ ] Rename `compare_foldx_evmutation.py` → `compare_spurs_eve.py`
- [ ] Re-derive quadrant thresholds — SPURS does not report kcal/mol the way FoldX did, so the old
      numeric cutoffs (ΔΔG > 1.0 destabilizing, etc.) do not transfer directly; pick new thresholds
      once real SPURS output exists (depends on T1.2 having run for real, which depends on Tier 2
      Phase 3)
- [ ] Update `integrate.smk`'s figure-generation target to reflect the new axis names

**Verification:** Run against synthetic/placeholder data first to confirm the merge logic (join on
family/position/ancestral_aa/derived_aa) still works; real quadrant classification waits on T1.2's
real run.
**Pass/fail:** Quadrant labels (`functional_gain`, `stability_function_tradeoff`, `deleterious`,
`neutral_drift`) are produced without crashing on a real convergent-sites TSV.

### T1.6 — Full-repo grep sweep for orphaned references
**Owner:** Claude Code | **Dependency:** T1.1–T1.5 all complete | **Effort:** ~1 hour
- [ ] `grep -rn "FoldX\|EVmutation\|EVcouplings" --include="*.md" --include="*.yaml" --include="*.py" --include="*.smk" .`
      outside `audit/`, `archive/`, and explicit historical-context sentences
- [ ] `grep -rn "esm2" config/config.yaml CLAUDE.md README.md` — confirm ESM-2 is only referenced
      as the Tier-3 baseline, never as primary

**Verification:** the grep commands above.
**Pass/fail:** Zero hits presenting the old stack as current/primary.

**Ledger:** - [ ] Tier 1 complete (method-stack migration)

---

## Tier 2 — Stub-phase implementation (build order matters — later phases depend on earlier ones)

### T2.1 — Structure prediction (`predict_af3.py`, `predict_chai1.py`, `assess_structure.py`, `classify_positions.py`)
**Owner:** Claude Code | **Dependency:** Tier 0, Tier 1 | **Effort:** ~2 days
- [ ] `predict_af3.py` — two-stage (CPU MSA → GPU inference), 5 seeds Tier 1 families
- [ ] `predict_chai1.py` — fallback
- [ ] `assess_structure.py` — pLDDT extraction, TM-align vs. reference PDB
- [ ] `classify_positions.py` — active site / pocket / surface / buried / interface / disordered
- [ ] Fill in `predict_structure.smk`, `structural_align.smk`

**Verification:** CLAUDE.md Test Gate 3 — AF3 vs. 7ZVA TM-score >0.90, Chai-1 fallback >0.85.
**Pass/fail:** as stated in the gate. This is also the unblock for T1.3 (ProSST/SaProt) — flag it
loudly if it fails, since two other tasks are waiting on it.

### T2.2 — Docking (`prepare_docking.py`, `run_docking.py`, `parse_docking.py`)
**Owner:** Claude Code | **Dependency:** T2.1 | **Effort:** ~3 days
**Verification:** CLAUDE.md Test Gate 5 (PQPQLPYP into neprosin 7ZVC, proline portion within 4 Å of
catalytic glutamates, visual + distance check).

### T2.3 — Electrostatics (`run_electrostatics.py`)
**Owner:** Claude Code | **Dependency:** T2.1 | **Effort:** ~2 days
**Verification:** pH 2.5/3.5/5.0 surface potential correlates directionally with species pitcher pH
for at least one family.

### T2.4 — Expression (`quantify_expression.py`)
**Owner:** Claude Code | **Dependency:** none (independent of structure work) | **Effort:** ~2 days
**Verification:** TPM matrix produced for the 3 named datasets (N. gracilis, Cephalotus, Dionaea).

### T2.5 — Integration (`build_atlas.py`, `generate_figures.py`, `map_convergence.py`)
**Owner:** Claude Code | **Dependency:** everything above + T1.5 | **Effort:** ~3 days
**Verification:** SQLite atlas builds without error; all 10 figures render, including the SPURS×
ProSST/SaProt/EVE quadrant plot (Fig. 5 equivalent) for at least one Tier-1 family.

### T2.6 — Webapp
**Owner:** Claude Code | **Dependency:** T2.5 | **Effort:** ~2 days
**Verification:** `streamlit run webapp/app.py` loads without error and can browse at least one
family's results.

**Ledger:** - [ ] Tier 2 complete

---

## Tier 3 — Manual/external data retrieval

Carried forward unchanged from the archived TASKS.md Tier 2 — not affected by the method-stack
pivot. See `archive/TASKS.md` for the detailed steps (BRAKER2 genome annotation, ENA/JGI manual
downloads, CDD domain verification) if this tier is picked up; re-verify accession numbers are
still current before running, since months may have passed.

- [ ] BRAKER2 on the 5 Tarnita 2023 genomes (Hellbender)
- [ ] ENA download of *N. mirabilis* proteins
- [ ] JGI download of *U. gibba* gene models
- [ ] CDD domain verification for accessions flagged "hypothetical protein"
- [ ] `local:` accession-prefix support in `fetch_sequences.py` (unblocks the three retrieval tasks
      above)

**Ledger:** - [ ] Tier 3 complete

---

## Tier 4 — Real test coverage

Replace the 5 placeholder test files, prioritized by whatever Tier 1/2 work lands first (a test
for a script that doesn't exist yet is wasted effort).

- [ ] `test_convergence.py` — real fixtures for `detect_convergence.py` (highest priority — this is
      the most-patched, most bug-prone script in the repo)
- [ ] `test_foldx_parse.py` → rename/rewrite as `test_spurs_parse.py` once T1.2 lands
- [ ] `test_evmutation.py` → rename/rewrite as `test_eve.py` once T1.4 lands
- [ ] `test_fetch.py` — real fixtures for `fetch_sequences.py`
- [ ] `test_docking_parse.py` — once T2.2 lands

**Ledger:** - [ ] Tier 4 complete

---

## Tier 5 — Manuscript

**Entry condition:** the SPURS×ProSST/SaProt/EVE quadrant figure (T2.5) has landed for at least one
Tier-1 family. No wet-lab section (Outcome A, PROJECT_PLAN.md §2.1).

- [ ] Draft section headings (IMRaD, MBE/GBE-tier)
- [ ] Cite the tightened novelty precedent set from PROJECT_PLAN.md §5 (Gerasimavicius/Livesey/
      Marsh, toxin-resistance-convergence papers) rather than CLAUDE.md's old overstated claim
- [ ] Distinguish CarnivorEnzyme from the nearest precedents: multi-family + cross-lineage +
      convergence-specific quadrant application (not the method itself, which is standard — see
      PROJECT_PLAN.md §5)

**Ledger:** - [ ] Tier 5 complete

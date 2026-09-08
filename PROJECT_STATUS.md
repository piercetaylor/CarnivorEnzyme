# CarnivorEnzyme — Project Status

> Date: 2026-08-26 (method-stack rows revised 2026-09-01 — see `audit/15_stability_predictor_audit.md`;
> first Hellbender run appended 2026-09-07 as §0 — see `audit/17_tier0_hellbender_verification.md`)
> Supersedes: `STATUS_2026-05-12.md` (archived, see `archive/STATUS_2026-05-12.md`)
> Companion documents: [PROJECT_PLAN.md](PROJECT_PLAN.md), [TODO.md](TODO.md)
> **This is a snapshot, not a living document.** It should be regenerated (not hand-patched into
> drift) at each major checkpoint — after a TODO.md tier closes, after an HPC run, after the
> method-stack migration lands. If you're reading this more than a few weeks after the date above,
> verify the phase/script table against the repo directly before trusting it.

---

## 0. First Hellbender run — 2026-09-07 (read this first)

The phylogeny→convergence chain ran on Hellbender for the first time. **It did not pass.** Full
detail in [audit/17_tier0_hellbender_verification.md](audit/17_tier0_hellbender_verification.md);
the parts that change how you should read the rest of this document:

- **T0.1 is still open.** 5 of 7 rules succeeded on real cluster software (MAFFT 7.525 → trimAl →
  IQ-TREE 3.1.3 WAG+G4/1000 UFBoot → root_tree, 18 sequences × 268 columns) and the chain then
  halted at `ancestral_reconstruction`. **`detect_convergence.py` never executed**, so the
  `9b1a8c5` node-matching fix — the thing Tier 0 exists to verify — is *still* untested against
  real data. §1's "top blocker" bullet remains accurate.
- **`config.yaml`'s `phylogeny.binary: iqtree3` is confirmed correct** and needs no change.
  Bioconda installs IQ-TREE 3.1.3 as `iqtree3` with `iqtree` as a symlink to it.
- **The halt is a paralogy finding, not a bug.** The four designated outgroup tips are not
  monophyletic in the gene tree — `Solanum_lycopersicum` sits sister to `Sarracenia_purpurea`,
  `Secale_cereale` sister to the two `Dionaea_muscipula` sequences — because the alignment holds 18
  sequences from only 9 species. `run_ancestral.py` refused to write a tree it could not orient,
  which is the guard working as designed. Agreed direction is to re-frame toward **recovering as
  many common ancestors as possible** (T0.1a–T0.1c in TODO.md) rather than forcing monophyly;
  family re-scoping to a true orthogroup stays separate, still-open work.
- **No pipeline logic was changed** by this run. New files are operational only: the SLURM driver
  `workflow/scripts/run_snakemake.sbatch`, the pass-criteria validator
  `workflow/scripts/validate_convergence_output.py`, `docs/HELLBENDER_RULES.md`, and `.claude/`
  guardrail hooks.
- **Running the pipeline on Hellbender needs setup the repo did not previously record.** The
  cluster's own Snakemake module cannot submit to SLURM, and conda envs cannot be built on the lab
  VAST share. See [docs/HELLBENDER_RULES.md](docs/HELLBENDER_RULES.md) before attempting a run.

### Note for the laptop: `CLAUDE.md` is now tracked

`CLAUDE.md` and `.claude/` were previously gitignored, so the laptop and the cluster each carried
their own untracked copy. As of 2026-09-07 the ignore rules are relaxed: `CLAUDE.md`,
`.claude/settings.json` and `.claude/hooks/` are **tracked**, and `results/sequences/` (50 FASTAs,
~58 KB) is tracked too so a clone is runnable without an rsync.

**Action needed on the laptop:** the cluster now has a `CLAUDE.md` scoped to *operating on
Hellbender*. The laptop's copy holds the science rationale and the Method Justification Table. The
intent is a **single unified, scoped `CLAUDE.md`** — merge the laptop's content in rather than
letting `git pull` overwrite it or spawning a second divergent file. Pull carefully the first time.
Note the Method Justification Table is stale on the FoldX / ESM-2 / EVmutation rows (§2), so the
merge is a chance to fix that.

---

## 1. Executive summary

CarnivorEnzyme is a computational atlas of convergently evolved digestive enzymes across
independently-evolved carnivorous plant lineages. As of this snapshot:

- **Phase 1 (retrieval), Phase 2 (convergence detection), Phase 3A (evolutionary-fitness triage),
  Phase 4A (ancestral structure), Phase 5D (FEP), Phase 5E (CpHMD)** have real implementations.
- **Phase 3 (structure prediction), Phase 4 (stability ΔΔG), Phase 5 (docking), Phase 6
  (electrostatics), Phase 7 (expression), Phase 8 (integration/atlas/figures)** are 100% stub
  scripts (`NotImplementedError` only) with matching stub `.smk` rule files.
- The method stack for stability, evolutionary-constraint, and DCA scoring changed on 2026-08-25
  (see [PROJECT_PLAN.md](PROJECT_PLAN.md) §2) — **FoldX, ESM-2-as-primary, and
  EVmutation/EVcouplings are no longer the target tools**, even though they're what's currently
  wired into `config.yaml` and CLAUDE.md's Method Justification Table. That migration is
  [TODO.md](TODO.md) Tier 1 and has **not yet been executed** — treat CLAUDE.md's tool table as
  stale on those three rows specifically until Tier 1 closes.
- The most recent architectural fix (ancestral-reconstruction node-matching, commit `9b1a8c5`,
  2026-08-24) has **not been verified against real Hellbender execution or real project data** —
  only synthetic data and a local Windows IQ-TREE binary. This is the current top blocker
  ([TODO.md](TODO.md) Tier 0).
- Test suite passes clean (13/13) but only 8 of those 13 tests exercise real logic
  (`test_config_consistency.py`); the other 5 test files are single-assertion placeholders.

---

## 2. Method stack (post-2026-08-25 decision, target state — not yet fully implemented)

| Axis | Tool | Status |
|------|------|--------|
| Structure prediction | AlphaFold3 (primary) / Chai-1 (fallback) / Boltz-2 (RNA-protein) | Config wired (`config.yaml structure:`); scripts are stubs |
| Stability ΔΔG | **FoldX 5.1** (Delgado et al. 2025) primary + **RaSP** (Blaabjerg et al. 2023) cross-check; SPURS (Li & Luo 2025) optional, not sole | **Not yet built** — `run_foldx_*.py` stubs are to be **completed, not deleted** (reversed 2026-09-01); `score_rasp.py` to be added. FoldX 5.1 is licensed and installed |
| Evolutionary fitness | **ProSST** (primary) + **SaProt** (secondary); ESM-2 retained as Tier-3 sanity baseline only | ESM-2 baseline implemented (`score_esm2.py`); ProSST/SaProt scripts don't exist yet, and need Foldseek 3Di tokens from a predicted structure — **hard-blocked on Phase 3 landing first** |
| Evolutionary coupling (DCA) | **EVE** (Frazer et al. 2021) | Stub scripts currently target EVmutation/EVcouplings (`run_evcouplings.py`, `run_evmutation.py`) — need renaming/rebuilding for EVE, not yet done |
| pH-dependent physics | **Electrostatics** (PDB2PQR + PROPKA + APBS) is the primary pH axis as of 2026-09-01; CpHMD is exploratory only | Electrostatics is a **stub** (`run_electrostatics.py`) and is now on the critical path. CpHMD is implemented but downgraded: it runs only on an unreleased GROMACS fork whose README asks users not to publish from it, and 3 pH points cannot yield a pKa |
| Binding free energy | Alchemical FEP (GROMACS + pmx) | Implemented (`run_fep.py`, `parse_fep.py`) |
| QM/MM | ORCA 6 + GROMACS | Implemented, gated off by default (`config.qmmm.enabled: false`) |

**FoldX is reinstated as the primary stability predictor** (reversed 2026-09-01; it had been slated
for removal on 2026-08-25). The reversal rests on three verified findings: FoldX is the one
predictor published as retaining correlation at *surface* residues, which is where Fukushima 2017
places the convergent substitutions; CpHMD computes protonated fractions and pKa, not folding free
energy, so it never was a substitute; and the benchmark figure used to justify the removal
("PCC ~0.40 on Megascale per MAVISp's March 2026 evaluation") does not exist in any source. The
stability axis is also re-scoped from thresholded kcal/mol to within-family rank, because direct
measurement on 7ZVA shows the thresholds are degenerate at surface sites. See
[audit/15_stability_predictor_audit.md](audit/15_stability_predictor_audit.md) and
[PROJECT_PLAN.md](PROJECT_PLAN.md) §2.2.

**Several citations in this repo are fabricated.** The 2025 FoldX revision is attributed throughout
CLAUDE.md and README.md to a nonexistent "Botte et al." (real: Delgado J, Reche R, Cianferoni D,
et al., *Bioinformatics* 41(2):btaf064); both GROMACS CpHMD references are invented, sharing one
bogus page number; and the archived May-2026 docs invent authors for SPURS and JanusDDG. In every
case the DOI is right and the author name is wrong — the signature of citations assembled from
search snippets. A Crossref sweep of all remaining citations is TODO Tier 1.

## 3. Phase-by-phase implementation table

Verified by direct inspection (line counts, not stub-template detection) on 2026-08-25.

| Phase | Scripts | Status | Matching `.smk` |
|-------|---------|--------|------------------|
| 1 — Retrieval | `fetch_sequences.py` (369 ln) | **Implemented** | `retrieve.smk` — implemented |
| 2 — Convergence | `detect_convergence.py` (742 ln); `root_tree.py` (214 ln), `run_ancestral.py` (317 ln) — both added by the Aug-24 ASR architecture fix | **Implemented** | `orthology.smk`, `alignment.smk`, `phylogeny.smk`, `convergence.smk` — all implemented |
| 3A — Evolutionary triage | `score_esm2.py` (269 ln) | **Implemented** (baseline only, post-decision) | `esm2.smk` — implemented; will be renamed `plm_scoring.smk` under TODO Tier 1 |
| 3 — Structure prediction | `predict_chai1.py`, `predict_af3.py`, `assess_structure.py`, `classify_positions.py` | **Stub** (`NotImplementedError` only) | `predict_structure.smk`, `structural_align.smk` — stub |
| 4A — Ancestral structure | `extract_ancestor.py` (269 ln), `compare_ancestor_modern.py` (218 ln) | **Implemented** — but `extract_ancestor.py` has the same node-naming-fallback pattern class just fixed in `detect_convergence.py`, not yet audited | `ancestral_structure.smk` — implemented |
| 4 — Stability ΔΔG | `run_foldx_repair.py`, `run_foldx_scan.py`, `parse_foldx.py` | **Stub — to be completed** (reversed 2026-09-01: no longer slated for removal). `score_rasp.py` to be added | `foldx.smk` — stub, to be filled in |
| DCA / evolutionary coupling | `run_evcouplings.py`, `run_evmutation.py`, `compare_foldx_evmutation.py` | **Stub — targets the wrong tool** (needs EVE rename/rebuild) | `evmutation.smk` — stub |
| 5 — Docking | `prepare_docking.py`, `run_docking.py`, `parse_docking.py` | **Stub** | `docking.smk` — stub |
| 5D — FEP | `run_fep.py` (294 ln), `parse_fep.py` (246 ln) | **Implemented** | `fep.smk` — implemented |
| 5E — CpHMD | `run_cphmd.py` (305 ln), `parse_cphmd.py` (262 ln) | **Implemented** | `cphmd.smk` — implemented |
| 5F — QM/MM | `run_qmmm.py` (286 ln) | **Implemented**, gated off by default | `qmmm.smk` — implemented |
| 6 — Electrostatics | `run_electrostatics.py` | **Stub** | `electrostatics.smk` — stub |
| 7 — Expression | `quantify_expression.py` | **Stub** | `expression.smk` — stub |
| 8 — Integration | `build_atlas.py`, `generate_figures.py`, `map_convergence.py` | **Stub** | `integrate.smk` — stub (covers both the SPURS×EVE quadrant comparison and the final atlas/figures) |
| Web interface | `webapp/app.py` (262 B) | **Stub** | `webapp/pages/*.py` — all 4 files are literally 0 bytes |

**Summary:** 9 of 22 named scripts implemented (plus 4 more added by the ASR fix: `root_tree.py`,
`run_ancestral.py`), 13 are stubs. 10 of 18 `.smk` rule files implemented, 8 are stubs. Config
sections exist for every phase, including the stub ones (`foldx:`, `docking:`, `electrostatics:`,
`expression:` in `config.yaml`) — config is consistently ahead of implementation.

---

## 4. Test suite status

```
pytest -q   →   13 passed
```

Don't read this as "well tested." Only `test_config_consistency.py` (8 assertions: accession/species
consistency, carnivory-origin counts, tier1/methods_benchmark disjointness, species-code
uniqueness, etc.) exercises real logic. `test_convergence.py`, `test_docking_parse.py`,
`test_evmutation.py`, `test_fetch.py`, `test_foldx_parse.py` each contain exactly one
`test_placeholder(): pass` — they pass trivially and cover nothing. Real coverage for the
implemented scripts (`detect_convergence.py`, `fetch_sequences.py`, `run_ancestral.py`, etc.) does
not exist yet. See [TODO.md](TODO.md) Tier 4.

---

## 5. Known open bugs and risks

| # | Issue | Blocks |
|---|-------|--------|
| 1 | Ancestral-reconstruction architecture fix (commit `9b1a8c5`, 2026-08-24) verified only on synthetic data + local Windows IQ-TREE binary — never run on Hellbender or real project alignments | Everything downstream of Phase 2 convergence detection — this is [TODO.md](TODO.md)'s Tier 0 blocker |
| 2 | Pre-existing Snakemake `ChildIOException`: `retrieve.smk`'s `fetch_all_sequences` rule outputs `results/sequences/`, which is the parent directory of `combine_family_sequences`'s output — Snakemake rejects a directory output that's an ancestor of another rule's output | Any real `snakemake -n` dry run past Phase 1 |
| 3 | `extract_ancestor.py` has the same class of node-naming fallback pattern (levelorder-index naming, substring `.state` matching) that was just fixed in `detect_convergence.py` — currently unreachable in the normal path post-fix, but not yet scrutinized | Trust in Phase 4A output before it's used for anything downstream |
| 4 | `config/enzyme_families.yaml.tmp.41128.1776193383614` — stray uncommitted temp file (21,855 B, dated Apr 14) sitting next to the real `enzyme_families.yaml` | Risk of confusion/accidental edit of the wrong file; no functional block |
| 5 | `hellbenderModules.txt` is 0 bytes | HPC module-load automation, if anything currently depends on it |
| 6 | `webapp/` is entirely stub — `app.py` is a 262-byte placeholder, all 4 `webapp/pages/*.py` files are literally empty | Phase 9 (web interface) has no runnable starting point |
| 7 | No `manuscript/` directory exists despite CLAUDE.md describing one | Manuscript work (Tier 5) has no scaffold yet |
| 8 | ~25 inline `TODO`/`NOT FOUND`/`PARTIAL_RESOLVE` sentinel accession entries remain in `config/enzyme_families.yaml` | Full-coverage retrieval for a handful of species/family combinations — already handled gracefully by `fetch_sequences.py`'s `_is_placeholder()` and caught by `test_config_consistency.py`, so not silently dangerous, just incomplete |

---

## 6. Documentation hygiene

- **Audit trail:** 15 numbered docs + `CHANGELOG.md` under `audit/`, covering the 2026-08-22
  (T1–T4) and 2026-08-24 (T5–T14) bug-fix/hardening/correction campaign. That work fixed silent
  zero-output failure modes in `detect_convergence.py`, `fetch_sequences.py`, and
  `04_run_hmmer_scan.sh`, plus taxonomic corrections (carnivory-origin count 3→1, MEROPS codes,
  GH19 Class I/IV split). See `audit/CHANGELOG.md` for the full index — not reproduced here.
- **This snapshot replaces** `STATUS_2026-05-12.md`, `TASKS.md`, `PROGRESS.md`,
  `temporary_TODO_42526.md`, and `PROJECT_RESTRUCTURE_2026-05-12.md` — all five described a method
  stack or task list that had drifted out of sync with the actual repo (either stale on the
  Aug-22–24 taxonomic corrections, or proposing an unexecuted pivot). Four are preserved under
  `archive/` with a superseded-by header. **`PROGRESS.md` is left at the repo root, untouched** —
  it's explicitly gitignored (`.gitignore:111`, the same convention CLAUDE.md itself uses) and was
  never tracked, so it's treated as scratch/local by repo convention; it's superseded in spirit by
  this document but wasn't moved into version control to respect that convention.
- **Second archival pass, 2026-09-01.** Six more root-level docs moved to `archive/` with
  superseded-by headers: `FOLDX_REVIEW_2026-05-12.md`, `FOLDX_ALTERNATIVES_AF3_2026-05-12.md`,
  `LITERATURE_REVIEW_2026-05-12.md`, `MANUSCRIPT_FRAMING_2026-05-12.md`,
  `ORTHOLOGY_AND_FAMILY_SCOPE_2026-05-12.md`, `ENZYME_EVOLUTION_PAPER_STRUCTURE.md`. The first two
  contain fabricated citations; the rest were stale on the Aug-22 origin-count correction or
  contradicted PROJECT_PLAN.md on the method stack. Root is now AGENTS.md, CLAUDE.md, PROGRESS.md
  (gitignored), PROJECT_PLAN.md, PROJECT_STATUS.md, README.md, TODO.md. **Two live findings were
  carried forward** out of `ENZYME_EVOLUTION_PAPER_STRUCTURE.md` rather than lost with it: the
  project has no dN/dS selection analysis and no PGLS phenotype correlation, both of which
  reviewers in this field expect (TODO.md T1.2d).
- **AGENTS.md is stale** — its status snapshot still says "As of April 2026… all 24 scripts are
  stubs", which PROJECT_STATUS.md §3 contradicts. Flagged, not rewritten.
- **README.md** still links to the old doc set and has its own known drift (e.g., a stale "nine
  independent lineage origins" reference) — out of scope for this documentation pass, but tracked
  as a follow-up in [TODO.md](TODO.md) Tier 0.5.

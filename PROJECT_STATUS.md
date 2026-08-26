# CarnivorEnzyme — Project Status

> Date: 2026-08-26
> Supersedes: `STATUS_2026-05-12.md` (archived, see `archive/STATUS_2026-05-12.md`)
> Companion documents: [PROJECT_PLAN.md](PROJECT_PLAN.md), [TODO.md](TODO.md)
> **This is a snapshot, not a living document.** It should be regenerated (not hand-patched into
> drift) at each major checkpoint — after a TODO.md tier closes, after an HPC run, after the
> method-stack migration lands. If you're reading this more than a few weeks after the date above,
> verify the phase/script table against the repo directly before trusting it.

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
| Stability ΔΔG | **SPURS** (Li & Luo, *Nat Commun* 17:891, 2025/2026) — AF-native, single axis | **Not yet built** — `score_spurs.py` doesn't exist yet; `run_foldx_*.py` stubs are slated for removal, not completion |
| Evolutionary fitness | **ProSST** (primary) + **SaProt** (secondary); ESM-2 retained as Tier-3 sanity baseline only | ESM-2 baseline implemented (`score_esm2.py`); ProSST/SaProt scripts don't exist yet, and need Foldseek 3Di tokens from a predicted structure — **hard-blocked on Phase 3 landing first** |
| Evolutionary coupling (DCA) | **EVE** (Frazer et al. 2021) | Stub scripts currently target EVmutation/EVcouplings (`run_evcouplings.py`, `run_evmutation.py`) — need renaming/rebuilding for EVE, not yet done |
| pH-dependent structural physics | **CpHMD** (GROMACS λ-dynamics via phbuilder) | Implemented (`run_cphmd.py`, `parse_cphmd.py`) — this is now the pipeline's *only* pH-resolved stability signal, since FoldX (the prior pH axis) is being dropped |
| Binding free energy | Alchemical FEP (GROMACS + pmx) | Implemented (`run_fep.py`, `parse_fep.py`) |
| QM/MM | ORCA 6 + GROMACS | Implemented, gated off by default (`config.qmmm.enabled: false`) |

**FoldX is being removed from the pipeline entirely**, including from what was CLAUDE.md's Test
Gate 4 (crystal-vs-prediction ΔΔG correlation check). See [PROJECT_PLAN.md](PROJECT_PLAN.md) §2
for the full rationale (SPURS is the newest peer-reviewed AF-native method; FoldX's pH-awareness
is real but redundant with CpHMD, which is already implemented and more rigorous; FoldX's own
general-ΔΔG benchmark is weak, PCC ~0.40–0.71 depending on dataset).

---

## 3. Phase-by-phase implementation table

Verified by direct inspection (line counts, not stub-template detection) on 2026-08-25.

| Phase | Scripts | Status | Matching `.smk` |
|-------|---------|--------|------------------|
| 1 — Retrieval | `fetch_sequences.py` (369 ln) | **Implemented** | `retrieve.smk` — implemented |
| 2 — Convergence | `detect_convergence.py` (742 ln); `root_tree.py` (214 ln), `run_ancestral.py` (317 ln) — both added by the Aug-24 ASR architecture fix | **Implemented** | `orthology.smk`, `alignment.smk`, `phylogeny.smk`, `convergence.smk` — all implemented |
| 3A — Evolutionary triage | `score_esm2.py` (269 ln) | **Implemented** (baseline only, post-decision) | `esm2.smk` — implemented; will be renamed `plm_scoring.smk` under TODO Tier 1 |
| 3 — Structure prediction | `predict_chai1.py`, `predict_af3.py`, `assess_structure.py`, `classify_positions.py` | **Stub** (`NotImplementedError` only) | `predict_structure.smk`, `structural_align.smk` — stub |
| 4A — Ancestral structure | `extract_ancestor.py` (269 ln), `compare_ancestor_modern.py` (218 ln) | **Implemented** — but `extract_ancestor.py` has the same node-naming-fallback pattern class just fixed in `detect_convergence.py`, not yet audited | `ancestral_structure.smk` — implemented |
| 4 — Stability ΔΔG | `run_foldx_repair.py`, `run_foldx_scan.py`, `parse_foldx.py` | **Stub — and slated for removal, not completion** (SPURS replaces this phase; see §2) | `foldx.smk` — stub, will not be filled in |
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

- **Audit trail:** 14 numbered docs + `CHANGELOG.md` under `audit/`, covering the 2026-08-22
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
- **README.md** still links to the old doc set and has its own known drift (e.g., a stale "nine
  independent lineage origins" reference) — out of scope for this documentation pass, but tracked
  as a follow-up in [TODO.md](TODO.md) Tier 0.5.

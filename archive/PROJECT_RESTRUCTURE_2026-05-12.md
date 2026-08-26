> **ARCHIVED 2026-08-25.** The open decisions in §8 of this document were resolved on 2026-08-25 —
> see `PROJECT_PLAN.md` §2 at the repo root for the resolutions and rationale. Preserved for
> history — do not treat the proposals here as current or as still-open.

# CarnivorEnzyme — Project Restructure

> Date: 2026-05-12
> Companion to [ENZYME_EVOLUTION_PAPER_STRUCTURE.md](ENZYME_EVOLUTION_PAPER_STRUCTURE.md), [STATUS_2026-05-12.md](STATUS_2026-05-12.md), [TASKS.md](TASKS.md)
> Purpose: Pivot the pipeline from method-organized phases to question-organized phases that mirror published enzyme evolution papers (Archetype B; Fukushima 2017 / Rubisco 2025 / Wang-Stiffler 2025 pattern). Define the specific publishable question.

---

## 1. The Specific Question (Real Outcome)

The pipeline now answers one specific question with a concrete testable hypothesis:

> **"Across nine independent origins of carnivory, do convergent amino acid substitutions in digestive enzymes track adaptation to pitcher fluid pH (1.9–6.0), and which evolutionary trajectory — functional gain, stability-function tradeoff, or neutral drift — characterizes each convergent site?"**

### The testable hypotheses

| H | Claim | Test | Pass criterion |
|---|-------|------|----------------|
| H1 | Convergent substitutions are biased toward acid-stabilizing positions | FoldX ΔΔG at pH 2.5 vs pH 7 across convergent sites vs random sites | Wilcoxon p < 0.01, effect size meaningful |
| H2 | Convergent substitution pattern correlates with species pitcher pH | PGLS regression of mean ΔΔG per species against pitcher pH | Lambda-controlled p < 0.05 |
| H3 | Convergent sites show positive selection on carnivorous lineages | PAML branch-site test M2a vs M1a; HyPhy aBSREL | LRT p < 0.05 at convergent sites; ω > 1 on carnivorous branches |
| H4 | Stability-function tradeoff is the dominant trajectory | Quadrant distribution: % sites in (FoldX-destab, evo-favored) quadrant vs other quadrants | > 40% in tradeoff quadrant for ≥ 3 of 5 families |
| H5 | (Optional, wet lab) The MRCA enzyme is inactive at pH 2.5 while modern is active | kcat/Km of resurrected ancestor vs modern at pH 2.5 vs 7 | Modern shows >10× kcat/Km advantage at pH 2.5 |

### Three landing scenarios

- **Outcome A (default — fully computational):** H1 + H2 + H3 + H4 supported. Target: *MBE* or *GBE*. Achievable 2026-08-15.
- **Outcome B (medium ambition):** A + H5 (one wet-lab figure). Target: *MBE* (lock) or *NEE* (probable). Achievable 2026-11-30.
- **Outcome C (high ambition):** B + cross-family trajectory pattern across 4–5 families. Target: *NEE* or *eLife*. Achievable 2027-Q1.

---

## 2. Restructured Phase Organization

The current Snakefile organizes phases by method (FoldX phase, FEP phase, etc.). The paper-organized structure reorganizes phases by **research-section purpose**:

### Old structure (method-driven)

```
Phase 1   Sequence retrieval
Phase 2   Alignment + phylogeny
Phase 2+  Orthology
Phase 3   Convergence detection
Phase 3A  ESM-2 triage
Phase 4   Structure prediction
Phase 4A  Ancestral structures
Phase 5A  FoldX
Phase 5B  EVcouplings
Phase 5C  FoldX × EVmutation
Phase 5D  Alchemical FEP            ← over-engineered
Phase 5E  Constant-pH MD            ← over-engineered
Phase 5F  QM/MM                     ← over-engineered (already gated)
Phase 6   Docking                   ← supplementary
Phase 7   Electrostatics            ← supplementary
Phase 8   Expression                ← supplementary
Phase 9   Atlas + 10 figures
```

### New structure (question-driven, paper-organized)

```
SECTION A — Comparative dataset
  A1. Sequence retrieval                      (retrieve.smk)
  A2. Genome annotation (Tarnita BRAKER2)     (braker2 scripts)

SECTION B — Phylogeny + Ancestral State
  B1. Alignment (MAFFT + trimAl)              (alignment.smk)
  B2. Phylogeny (IQ-TREE 3 + ASR)             (phylogeny.smk)
  B3. Tree rooting                            (phylogeny.smk → root_tree)
  B4. Ancestral structure prediction          (ancestral_structure.smk)

SECTION C — Selection + Convergence
  C1. Site-level convergence (Fukushima 2017)  (convergence.smk)  ✅
  C2. ωC error-corrected convergence           NEW: oc_convergence.smk
  C3. PLM embedding convergence (ACEP)         NEW: acep.smk
  C4. dN/dS site-level (PAML M1a/M2a)          NEW: selection.smk ← MISSING
  C5. dN/dS branch-site (PAML A; HyPhy)        NEW: selection.smk ← MISSING

SECTION D — Structure Prediction
  D1. AF3 / Chai-1                            (predict_structure.smk)
  D2. Structure assessment                     (assess_structure.py)
  D3. Active-site classification               (classify_positions.py)
  D4. Structural alignment (Foldseek)          (structural_align.smk)

SECTION E — Mutation Effect Prediction
  E1. FoldX ΔΔG at pH 2.5 / 3.5 / 5.0           (foldx.smk)
  E2. ML ΔΔG (SPURS or RaSP)                   NEW: stability_ml.smk
  E3. EVcouplings / EVmutation                  (evmutation.smk)
  E4. ESM-2 / ProSST LLR                        (esm2.smk → plm_scoring.smk)
  E5. Stability × Evolution quadrant            (integrate.smk → compare_stability_evolution.py)

SECTION F — Phenotype Correlation (NEW)
  F1. Pitcher pH dataset assembly              NEW: phenotype.smk
  F2. PGLS regression                          NEW: pgls.py
  F3. Continuous trait mapping                  NEW: trait_mapping.py

SECTION G — Integration
  G1. Atlas SQLite                              (integrate.smk → build_atlas.py)
  G2. Publication figures (8)                   (integrate.smk → generate_figures.py)

SECTION S — Supplementary / Preliminary Data (gated, runs only on top hits)
  S1. Constant-pH MD                            (cphmd.smk) — gate on E5 quadrant
  S2. Alchemical FEP                            (fep.smk) — ≤ 20 sites total
  S3. QM/MM                                     (qmmm.smk) — disabled by default
  S4. Substrate docking                         (docking.smk) — Vina ranking only
  S5. Electrostatics (APBS)                     (electrostatics.smk) — figure support
  S6. Expression (Salmon)                       (expression.smk) — figure support

SECTION W — Wet Lab (optional, gates Outcome B/C)
  W1. Ancestor + modern variant gene synthesis  external
  W2. Expression in Pichia / E. coli            external
  W3. Activity assay at pH 2.5 / 5 / 7          external
```

### Phase target convenience rules

| Target | Equivalent to | What it produces |
|--------|---------------|-------------------|
| `phase_dataset` | A complete | All sequence FASTAs in results/sequences/ |
| `phase_phylogeny` | B complete | Trees, .state files, ancestors |
| `phase_evolution` | C complete | Convergence + dN/dS TSVs |
| `phase_structure` | D complete | AF3 predictions + classification |
| `phase_effects` | E complete | FoldX + ML ΔΔG + LLR + quadrant |
| `phase_phenotype` | F complete | PGLS regression results |
| `phase_atlas` | G complete | SQLite + 8 figures |
| `phase_supplementary` | S complete | CpHMD/FEP/docking/electrostatics/expression |
| `manuscript` | All of A–G | The publishable computational package |
| `manuscript_plus` | manuscript + W | Outcome B (wet-lab augmented) |

---

## 3. Codebase Audit — Script-by-Script Mapping

### What's already there (15 implemented scripts)

| Script | Section | Status | Action |
|--------|---------|--------|--------|
| fetch_sequences.py | A1 | ✅ implemented | Keep; add `local:` prefix support (T1.5) |
| root_tree.py | B3 | ✅ implemented | Keep |
| detect_convergence.py | C1 | ✅ implemented | Keep |
| extract_ancestor.py | B4 | ✅ implemented | Keep |
| compare_ancestor_modern.py | B4 | ✅ implemented | Keep — already above-standard |
| score_esm2.py | E4 | ✅ implemented | Rename → `score_plm.py`, support ProSST as primary |
| run_cphmd.py | S1 | ✅ implemented | Move to Supplementary, gate strictly |
| parse_cphmd.py | S1 | ✅ implemented | Move to Supplementary |
| run_fep.py | S2 | ✅ implemented | Move to Supplementary, ≤ 20 sites cap |
| parse_fep.py | S2 | ✅ implemented | Move to Supplementary |
| run_qmmm.py | S3 | ✅ implemented (gated off) | Leave gated off |
| braker2/*.sh | A2 | ✅ implemented | Run on Hellbender (T2.1) |

### What's stubbed (15 stub scripts that need implementing)

| Script | Section | Priority | New script needed? |
|--------|---------|----------|---------------------|
| predict_af3.py | D1 | HIGH | Keep stub; implement |
| predict_chai1.py | D1 | MEDIUM | Keep stub; implement as fallback |
| assess_structure.py | D2 | HIGH | Keep stub; implement |
| classify_positions.py | D3 | HIGH | Keep stub; implement |
| run_foldx_repair.py | E1 | HIGH | Keep stub; implement |
| run_foldx_scan.py | E1 | HIGH | Keep stub; **add PypKa upstream** |
| parse_foldx.py | E1 | HIGH | Keep stub; implement |
| run_evcouplings.py | E3 | MEDIUM | Keep stub; implement when Neff ≥ 500 |
| run_evmutation.py | E3 | MEDIUM | Keep stub; implement |
| compare_foldx_evmutation.py | E5 | HIGH | Rename → `compare_stability_evolution.py`; implement 2D quadrant |
| prepare_docking.py | S4 | LOW | Keep stub; move to supplementary tier |
| run_docking.py | S4 | LOW | Keep stub; move to supplementary tier |
| parse_docking.py | S4 | LOW | Keep stub; move to supplementary tier |
| run_electrostatics.py | S5 | LOW | Keep stub; move to supplementary tier |
| quantify_expression.py | S6 | LOW | Keep stub; move to supplementary tier |
| build_atlas.py | G1 | HIGH | Keep stub; implement after E5 |
| generate_figures.py | G2 | HIGH | Keep stub; implement after F2 |
| map_convergence.py | D3 | LOW | Likely redundant with classify_positions.py — **remove** |

### What's missing — new scripts to create

| Script | Section | Priority | Notes |
|--------|---------|----------|-------|
| **`run_codeml.py`** | C4 | **HIGH** | PAML codeml M1a/M2a site model + LRT |
| **`run_absrel.py`** | C5 | **HIGH** | HyPhy aBSREL branch-site selection |
| **`run_busted.py`** | C5 | MEDIUM | HyPhy BUSTED whole-gene selection |
| **`run_relax.py`** | C5 | MEDIUM | HyPhy RELAX selection-relaxation test |
| **`pgls.py`** | F2 | **HIGH** | Phylogenetic generalized least squares (uses R `caper` or Python `phylolm`) |
| `run_oc.py` | C2 | LOW | ωC error-corrected convergence (Fukushima 2023) |
| `run_acep.py` | C3 | LOW | ACEP PLM embedding convergence (Wang 2025) |
| `score_spurs.py` | E2 | MEDIUM | AF-native ML ΔΔG (Foldetta-style second axis) |
| `score_rasp.py` | E2 | LOW | Alternative second-axis ΔΔG |
| `phenotype_assemble.py` | F1 | HIGH | Compile pitcher_pH per species from species.yaml; output trait table |
| `trait_mapping.py` | F3 | LOW | Visualize pitcher_pH on phylogeny (phytools or ete3) |

**Estimated work:** 5 high-priority new scripts (codeml, absrel, pgls, phenotype_assemble, score_spurs). Each ~half-day to 1 day. Total: ~5 days.

### What's over-engineered — move to Supplementary tier

| Currently | New tier | Justification |
|-----------|----------|---------------|
| Phase 5D FEP (run_fep.py + fep.smk) | Section S2 | No published enzyme evolution paper uses FEP; gate strictly on E5 quadrant winners; ≤ 20 sites total |
| Phase 5E CpHMD (run_cphmd.py + cphmd.smk) | Section S1 | No published enzyme evolution paper uses CpHMD; supports mechanism not phenotype; gate on E5 |
| Phase 5F QM/MM (run_qmmm.py + qmmm.smk) | Section S3 | Already gated off; keep gated off for v1 manuscript |
| Phase 6 docking | Section S4 | Butts 2016 precedent exists but secondary to pH adaptation question |
| Phase 7 electrostatics | Section S5 | Figure-only support; not a primary claim |
| Phase 8 expression | Section S6 | Only 3 datasets; tangential to convergence question |

---

## 4. Restructured Manuscript Figure Set

Old plan: 10 figures organized by method.

New plan: 8 main figures + supplementary, organized by paper section:

| Fig | Section | Content | Required for |
|-----|---------|---------|--------------|
| 1 | B + C | Per-family phylogeny with convergent sites highlighted on lineage branches | H1 framing |
| 2 | D | Structural superposition per family; convergent positions colored by active-site/surface/buried | H1 framing |
| 3 | E1 | FoldX ΔΔG distributions at pH 2.5, 3.5, 5.0 — convergent vs background | H1 test |
| 4 | E3 + E4 | EVmutation + ESM-2/ProSST LLR distributions for convergent sites | H1 + H4 |
| 5 | E5 | **FoldX × evolutionary quadrant per family — CORE NOVEL FIGURE** | H4 test |
| 6 | F2 | PGLS regression: mean ΔΔG_pH2.5 per species vs species pitcher_pH | **H2 test** |
| 7 | C4 + C5 | dN/dS site map + branch-site selection on carnivorous lineages | **H3 test** |
| 8 | (Optional) W | kcat/Km ancestor vs modern at pH 2.5 / 5 / 7 — single neprosin | H5 test (Outcome B/C only) |

### Supplementary figures (preliminary mechanism data)

| Sup Fig | Content |
|---------|---------|
| S1 | CpHMD pKa shifts at active-site Asp/Glu (top quadrant winners only) |
| S2 | FEP ΔΔG_bind for ≤ 20 priority sites |
| S3 | Substrate docking affinities |
| S4 | Electrostatic surfaces at pH 2.5 / 5 |
| S5 | Expression heatmap (3 datasets) |
| S6 | Ancestral vs modern structural displacement (compare_ancestor_modern output) |
| S7 | Per-family Neff and MSA quality |
| S8 | Convergence detection sensitivity analysis |

---

## 5. Code Reorganization Changes

### Files to rename

| Old | New | Why |
|-----|-----|-----|
| `workflow/scripts/score_esm2.py` | `workflow/scripts/score_plm.py` | Now supports ProSST as primary, ESM-2 as baseline |
| `workflow/scripts/compare_foldx_evmutation.py` | `workflow/scripts/compare_stability_evolution.py` | Generic across stability + evolution axes |
| `workflow/rules/esm2.smk` | `workflow/rules/plm_scoring.smk` | Consistent with renamed script |

### Files to delete (redundant)

| File | Why |
|------|-----|
| `workflow/scripts/map_convergence.py` | Functionality covered by classify_positions.py |

### New files to create

```
workflow/scripts/
├── run_codeml.py            # PAML site-level dN/dS
├── run_absrel.py            # HyPhy branch-site selection
├── pgls.py                  # Phylogenetic generalized least squares
├── phenotype_assemble.py    # Build pitcher_pH trait table from species.yaml
├── score_spurs.py           # AF-native ML ΔΔG (optional second axis)

workflow/rules/
├── selection.smk            # dN/dS pipeline rules
├── phenotype.smk            # PGLS pipeline rules
├── stability_ml.smk         # SPURS / RaSP (optional)

workflow/envs/
├── paml.yaml                # PAML codeml + ete3
├── hyphy.yaml               # HyPhy binaries
├── pgls.yaml                # R + caper + phylolm; or Python phylolm-equivalent
```

### Snakefile reorganization

The new Snakefile structure:

```python
configfile: "config/config.yaml"

# ── Rule modules organized by paper section ──────────────────────────────────
# Section A — Comparative dataset
include: "workflow/rules/retrieve.smk"
include: "workflow/rules/orthology.smk"

# Section B — Phylogeny + ASR
include: "workflow/rules/alignment.smk"
include: "workflow/rules/phylogeny.smk"
include: "workflow/rules/ancestral_structure.smk"

# Section C — Selection + Convergence
include: "workflow/rules/convergence.smk"
include: "workflow/rules/selection.smk"           # NEW: PAML + HyPhy

# Section D — Structure
include: "workflow/rules/predict_structure.smk"
include: "workflow/rules/structural_align.smk"

# Section E — Mutation Effects
include: "workflow/rules/foldx.smk"
include: "workflow/rules/stability_ml.smk"        # NEW: SPURS / RaSP
include: "workflow/rules/evmutation.smk"
include: "workflow/rules/plm_scoring.smk"         # renamed from esm2.smk

# Section F — Phenotype Correlation
include: "workflow/rules/phenotype.smk"           # NEW: PGLS

# Section G — Integration
include: "workflow/rules/integrate.smk"

# Section S — Supplementary (gated, downstream of E5)
include: "workflow/rules/cphmd.smk"
include: "workflow/rules/fep.smk"
include: "workflow/rules/qmmm.smk"
include: "workflow/rules/docking.smk"
include: "workflow/rules/electrostatics.smk"
include: "workflow/rules/expression.smk"

# ── Section-level convenience targets ────────────────────────────────────────
rule phase_dataset:
    input: "results/sequences/"

rule phase_phylogeny:
    input:
        expand("results/phylogenies/{family}.rooted.treefile", family=FAMILIES),
        expand("results/ancestral/{family}.mrca_ancestor.fa", family=FAMILIES),

rule phase_evolution:
    input:
        expand("results/convergence/{family}.convergent_sites.tsv", family=FAMILIES),
        expand("results/selection/{family}.codeml_results.tsv", family=FAMILIES),
        expand("results/selection/{family}.absrel_results.tsv", family=FAMILIES),

rule phase_structure:
    input:
        expand("results/structures/{family}.af3.cif", family=FAMILIES),
        expand("results/structures/{family}.classified_positions.tsv", family=FAMILIES),

rule phase_effects:
    input:
        expand("results/foldx/{family}.ddg_matrix.tsv", family=FAMILIES),
        expand("results/stability_ml/{family}.spurs_scores.tsv", family=FAMILIES),
        expand("results/evmutation/{family}.evmut_scores.tsv", family=FAMILIES),
        expand("results/plm/{family}.prosst_scores.tsv", family=FAMILIES),
        expand("results/integration/{family}.quadrant.tsv", family=FAMILIES),

rule phase_phenotype:
    input:
        "results/phenotype/pgls_results.tsv",
        "results/phenotype/trait_mapping.pdf",

rule phase_atlas:
    input:
        "results/atlas/atlas.sqlite",
        expand("results/atlas/figures/fig{n}.pdf", n=range(1, 9)),  # 8 main figures

rule phase_supplementary:
    """Optional: runs CpHMD/FEP/docking/electrostatics/expression on top hits only."""
    input:
        # gated outputs only
        expand("results/cphmd/{family}.pka_shifts.tsv", family=FAMILIES),
        expand("results/fep/{family}.fep_results.tsv", family=FAMILIES),
        expand("results/docking/{family}.top_affinities.tsv", family=FAMILIES),
        expand("results/electrostatics/{family}.surface_potentials.tsv", family=FAMILIES),
        expand("results/expression/{species}.target_tpm.tsv", species=EXPRESSION_SPECIES),

rule manuscript:
    """The publishable computational package (Outcome A)."""
    input:
        rules.phase_atlas.input,

rule manuscript_plus:
    """Outcome B — manuscript + supplementary mechanism data."""
    input:
        rules.manuscript.input,
        rules.phase_supplementary.input,
```

---

## 6. Migration Plan

### Step 1 — Documentation (this week)

- [x] Create [ENZYME_EVOLUTION_PAPER_STRUCTURE.md](ENZYME_EVOLUTION_PAPER_STRUCTURE.md)
- [x] Create this restructure document
- [ ] Update [README.md](README.md) with the specific question
- [ ] Update [STATUS_2026-05-12.md](STATUS_2026-05-12.md) §2 (Scientific Framing) to use the specific hypotheses
- [ ] Update [TASKS.md](TASKS.md) Tier 1: add selection.smk + pgls.py tasks; demote Supplementary section tasks

### Step 2 — Add the missing methods (1–2 weeks)

- [ ] Implement `run_codeml.py` and `selection.smk` for PAML site-level dN/dS (M1a vs M2a)
- [ ] Implement `run_absrel.py` for HyPhy branch-site selection on carnivorous lineages
- [ ] Implement `phenotype_assemble.py` to build pitcher_pH trait table from species.yaml
- [ ] Implement `pgls.py` for the PGLS regression of mean ΔΔG vs pitcher_pH

### Step 3 — Restructure the Snakefile (1 day)

- [ ] Rename rule files where needed
- [ ] Reorganize Snakefile to section-based includes
- [ ] Replace phase1–phase9 targets with phase_dataset / phase_phylogeny / phase_evolution / etc.
- [ ] Add `manuscript` and `manuscript_plus` final targets
- [ ] Move supplementary phases (CpHMD, FEP, docking, electrostatics, expression) into a single `phase_supplementary` rule with strict gating

### Step 4 — Implement the critical-path stubs (3–4 weeks)

Critical path: predict_af3.py → assess_structure.py → classify_positions.py → run_foldx_repair.py → run_foldx_scan.py → parse_foldx.py → score_plm.py (rename + ProSST upgrade) → compare_stability_evolution.py → pgls.py → build_atlas.py → generate_figures.py

### Step 5 — Run and analyze (4–6 weeks on Hellbender)

- AF3 prediction (~200 sequences, ~5 days HPC)
- FoldX at 3 pH values (~3 days)
- ML ΔΔG (SPURS, hours)
- EVmutation for deep-MSA families (~2 days)
- PLM scoring (~1 day)
- Quadrant analysis (hours)
- PGLS regression (hours)
- Figure generation (~1 day)

### Step 6 — Manuscript decision point

After phase_atlas completes:

- If H1 + H2 + H3 + H4 all hold: write Outcome A manuscript → bioRxiv + MBE submission
- If results are noisy or mixed: target plant journal (Plant Cell, New Phytologist)
- If results are strong: initiate wet-lab collaboration for Outcome B

### Step 7 — Optional Section S (Supplementary mechanism)

Only run S1–S6 after phase_effects completes and quadrant winners are identified:

- CpHMD on top 5 active-site Glu/Asp shifts per family
- FEP on ≤ 20 priority sites total
- Docking on top quadrant orthologs
- Electrostatics on representative species per origin
- Expression on the 3 available datasets

These appear in the manuscript as supplementary figures only.

---

## 7. What This Restructure Achieves

| Before | After |
|--------|-------|
| Atlas / resource paper (vague question) | Hypothesis-testing paper (5 explicit H1–H5) |
| 9 phases, method-organized | 7 paper-section phases + 1 supplementary |
| 10 figures, scattered | 8 main figures matching IMRaD structure |
| No dN/dS, no PGLS — reviewers will demand | dN/dS + PGLS as Section C and Section F |
| CpHMD/FEP/QM-MM as main pipeline (over-engineered) | Same tools as Supplementary mechanism data |
| Implicit target: NEE/Science | Explicit targets: MBE default, NEE with wet lab |
| Compute budget ~25k GPU-hours | Compute budget ~5–10k GPU-hours (supplementary gated) |

The restructure does not lose anything implemented. It re-tiers what's already built, adds the two missing standard methods (dN/dS + PGLS), and explicitly defines what the paper answers.

---

## 8. Open Decisions Required from PI

These need a decision before Step 2 begins:

1. **Pick the target journal tier.** Outcome A (MBE/GBE default) or commit to Outcome B (initiate wet-lab now)?
2. **Pick the second stability axis: SPURS vs RaSP.** SPURS is AF-native (best benchmark on AF inputs); RaSP is faster and gives saturation maps. Recommend SPURS.
3. **Pick the primary PLM: ProSST vs SaProt vs ESM-2.** Recommend ProSST (ProteinGym top); keep ESM-2 as baseline (already implemented).
4. **EVcouplings or EVE?** Same MSA requirements; EVE has better benchmark numbers but EVmutation has more published precedent. Recommend EVE.
5. **Wet-lab collaboration timeline.** If pursuing Outcome B, start outreach now (Goh lab Malaysia or Brodelius lab Sweden).

---

## 9. Citation Targets

The restructured project should cite, in Methods and Discussion:

- **Fukushima 2017** — sequence-level convergence framing
- **Fukushima & Pollock 2023** — convergence statistical method
- **Albert/Fukushima 2026** — open question framing
- **Pavlovič 2025** + **Freund 2022** — pitcher pH phenotype
- **Iñiguez 2022 Sci Adv (Solanaceae Rubisco)** — closest precedent for "ancestral reconstruction + plant enzyme + computational" workflow
- **Rubisco convergence bioRxiv 2025.10.08.681247** — closest computational precedent
- **Wang/Stiffler 2025 PNAS (ACEP)** — convergence detection complement
- **Gerasimavicius, Livesey, Marsh 2023** — Foldetta-style two-stability-method consensus
- **Botte 2025 Bioinformatics 41:btaf064** — FoldX 5.1 pH-aware force field
- **Yang 2007** — PAML codeml
- **Smith et al. 2015** — HyPhy aBSREL
- **Revell 2024** — phytools R package
- **Akdel 2022 NSMB** — FoldX-on-AlphaFold validation
- **Hopf 2017** — EVmutation
- **Lin 2023 Science (ESM-2)** OR **Li 2024 NeurIPS (ProSST)**

This citation profile matches both Archetype A (Iñiguez, Thornton-style) and Archetype B (Fukushima, Rubisco-style) — making the paper defensible from both angles.

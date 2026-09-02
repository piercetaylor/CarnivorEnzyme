> **ARCHIVED 2026-09-01.** Superseded by `audit/15_stability_predictor_audit.md`.
> Stale on the origin count ("nine") and contradicts PROJECT_PLAN.md §2 on CpHMD (it calls CpHMD/FEP/QM-MM over-engineered; PROJECT_PLAN made CpHMD the sole pH axis). Its two live findings — missing dN/dS selection tests and missing PGLS — are carried forward into the 2026-09-01 audit and TODO.md.
> Preserved for history — do not treat as current.

# Enzyme Evolution Paper Structure — Reference Template

> Date: 2026-05-12
> Purpose: Document the methods structure of typical enzyme evolution papers (2023–2026), framed around the open questions Fukushima et al. proposed. Use this template to audit and restructure CarnivorEnzyme to match how the field actually publishes.

---

## 1. The Fukushima Open Questions

The CarnivorEnzyme project was justified by open questions from three Fukushima-group publications. Re-stating them precisely:

### Q1 — Fukushima et al. 2017 (*Nat. Ecol. Evol.* 1:0059)

> "Convergent amino acid substitutions in digestive enzymes of independently-evolved carnivorous plant lineages occur at a rate exceeding neutral expectation. **Whether these substitutions confer specific structural, functional, or biochemical changes — and whether those changes are adaptive — remains to be tested.**"

The 2017 paper identified convergence at sequence level only. It did **not** test functional consequence, did **not** predict structures, did **not** run any biophysical analysis.

### Q2 — Fukushima & Pollock 2023 (*Nat. Ecol. Evol.* 7:1156)

> "Detection of convergent substitutions requires error correction (the ωC metric) because branch-length variation inflates false positives. **Whether ωC-significant convergent sites in carnivorous plant enzymes are mechanistically responsible for the digestive phenotype is an open question.**"

This paper provided the statistical framework. It did **not** apply it to functional analysis.

### Q3 — Albert, Fukushima et al. 2026 (*Trends in Genetics* doi:S0168-9525(26)00034-X)

The 2026 review explicitly names three open questions:

- **Q3a:** Functional consequences of convergent substitutions at the structural level
- **Q3b:** Polyploidy × convergence interaction — which subgenome copies retain convergent substitutions in decaploid *N. gracilis*, tetraploid *Dionaea*, dodecaploid *D. capensis*?
- **Q3c:** Tradeoffs between local and global gene duplication and intergenic DNA deletion

Q3a is the question CarnivorEnzyme is built to answer. Q3b is the polyploidy angle currently under-developed.

### Q4 — Pavlovič 2025 (*New Phytologist* doi:10.1111/nph.70229) + Freund 2022 (*Plant Physiology* 190:44)

Pitcher fluid pH spans pH 1.9 (*N. rafflesiana*) to pH 6 (*N. ampullaria*) across species and lineages. **Whether digestive enzyme sequence/structure tracks this pH variation is an open phylogenetic-comparative question.**

This is the pH-phenotype axis that gives CarnivorEnzyme a quantitative trait to correlate against.

---

## 2. The Two Archetypes of Published Enzyme Evolution Papers

The enzyme evolution literature follows two main methodological archetypes. Hybrid papers exist but are rare. Top-tier impact (Nature, Science, eLife) effectively requires Archetype A; computational-only papers cluster in Archetype B venues (MBE, GBE, Plant Cell, New Phytologist).

### Archetype A — "Resurrect and Assay"

**Method backbone:** Phylogeny → ASR → gene synthesis → expression → kinetic characterization → (optional) crystal structure.

**Canonical examples:**

- **Iñiguez et al. 2022** (*Sci. Adv.* doi:10.1126/sciadv.abm6871) — Resurrected ancestral Solanaceae Rubisco. Measured Vc, Sc/o in vivo. No FoldX.
- **McKeown & Thornton 2024** (*Nature*) — DMS on resurrected steroid receptor ancestors.
- **Anderson et al. 2020** (*PNAS*) — Hormone receptor ASR + ligand-binding kinetics.
- **Hochberg lab** — Livada et al. 2023 (*ACS Catal.*) ene reductase ASR + Tm.
- **JACS Au 2024** (doi:10.1021/jacsau.4c00653) — Field-defining review: expression + kinetic characterization is the default endpoint.

**Required figures:** Phylogeny, ASR posterior probabilities, Km/kcat tables, Tm curves, sometimes a crystal structure.

**Journal landing:** *Nature*, *Science*, *eLife*, *Nat. Ecol. Evol.*, *Nat. Chem. Biol.*

### Archetype B — "Compute and Correlate"

**Method backbone:** Phylogeny → ASR → convergence/selection test → structure prediction → ΔΔG mapping → phylogenetic comparative test against phenotype.

**Canonical examples:**

- **Fukushima 2017** (*Nat. Ecol. Evol.*) — Phylogeny + convergence counting. No structure, no ΔΔG.
- **Fukushima & Pollock 2023** (*Nat. Ecol. Evol.*) — ωC metric. Methods paper.
- **Rubisco convergence bioRxiv 2025.10.08.681247** — Phylogeny + ASR + FoldX on AF/spinach crystal + PGLS vs. growing-season temperature. Single family. **FoldX-only.** Published anyway.
- **Wang/Stiffler et al. 2025** (*PNAS* doi:10.1073/pnas.2418254122) — ACEP: PLM embeddings for adaptive convergence detection. No structure, no ΔΔG.
- **Butts et al. 2016** (*Proteins* 84:1517) — Homology modeling + docking, single lineage. Pre-AlphaFold.
- **Riziotis et al. 2025** (*FEBS J.* doi:10.1111/febs.17332) — 34-case review of convergent enzyme evolution. Methods used: sequence comparison, CATH fold matching, M-CSA active-site geometry, catalytic mechanism similarity. **No FoldX, FEP, CpHMD anywhere.**

**Required figures:** Phylogeny with convergent sites highlighted, structural mapping of convergent positions, ΔΔG or selection-test heatmap, phylogenetic comparative regression (PGLS) against the phenotype axis.

**Journal landing:** *Molecular Biology and Evolution* (MBE), *Genome Biology and Evolution* (GBE), *Plant Cell*, *Plant Physiology*, *New Phytologist*.

---

## 3. The Standard Methods Stack (Archetype B in detail)

This is the structure CarnivorEnzyme should mirror — minus over-engineering, plus the missing pieces.

### 3.1 Comparative Dataset Assembly

| Step | Standard tool | Notes |
|------|---------------|-------|
| Sequence retrieval | NCBI Entrez, UniProt API, ENA | Curate accessions; verify domain composition |
| Genome annotation (un-annotated) | BRAKER2 / AUGUSTUS | When only raw assembly exists |
| Outgroup selection | Phylogenetic neighbors of focal lineages | Critical for ancestral state |
| Manual curation | Per-accession verification | Filter wrong-family / wrong-organism entries |

### 3.2 Phylogeny

| Step | Standard tool | Notes |
|------|---------------|-------|
| Alignment | MAFFT L-INS-i | For < 200 sequences |
| Trimming | trimAl `-automated1` | Inspect HTML output |
| Model selection | IQ-TREE ModelFinder | Or PartitionFinder |
| ML inference | IQ-TREE with UFBoot2 + SH-aLRT | UFBoot ≥95 + SH-aLRT ≥80 = strong support |
| Rooting | Outgroup or midpoint | Required for ancestral state direction |

### 3.3 Ancestral Sequence Reconstruction

| Step | Standard tool | Notes |
|------|---------------|-------|
| Marginal ASR | IQ-TREE `-ancestral` → .state | Posterior probability per site per node |
| Joint ASR (alternative) | PAML codeml `-baseml` | Slower, sometimes preferred for resurrection |
| Posterior threshold | ≥ 0.8 | Sites below: mask as 'X' or flag uncertainty |
| MRCA identification | Manual or ete3 traversal | Tag the carnivore MRCA node per family |

### 3.4 Selection Analysis (**MISSING FROM CARNIVORENZYME — ADD**)

This is the single most expected element by reviewers that the project currently lacks.

| Step | Standard tool | Notes |
|------|---------------|-------|
| Site-level dN/dS | PAML codeml M1a/M2a/M7/M8 | Likelihood ratio test for positively-selected sites |
| Branch-site dN/dS | PAML codeml branch-site model A | Tests selection on specific branches (carnivorous lineages) |
| Episodic selection | HyPhy aBSREL | Branch-by-branch episodic selection |
| Pervasive selection | HyPhy BUSTED | Whole-gene positive selection |
| Selection relaxation | HyPhy RELAX | Tests for relaxed vs intensified selection |
| Convergence (site coincidence) | Fukushima 2017 method | Already implemented in detect_convergence.py |
| Convergence (error-corrected) | ωC Fukushima & Pollock 2023 | Cited but not yet implemented |
| Convergence (PLM embedding) | ACEP Wang 2025 | Adds an orthogonal axis |

### 3.5 Structure Prediction

| Step | Standard tool | Notes |
|------|---------------|-------|
| Primary predictor | AlphaFold3 (5 seeds Tier 1) | Or Chai-1 fallback |
| Quality assessment | pLDDT, PAE | Field standard: pLDDT ≥ 70 at sites of interest |
| Reference benchmarking | TM-align vs. crystal structure | TM-score > 0.90 target |
| Disorder masking | pLDDT < 50 excluded | Some labs use < 70 |

### 3.6 Mutation Effect Prediction

| Step | Standard tool | Notes |
|------|---------------|-------|
| Stability ΔΔG (physics) | FoldX 5.1 RepairPDB + PositionScan | At pH 2.5/3.5/5.0 if pH-relevant |
| Stability ΔΔG (ML, optional) | RaSP or SPURS | "Foldetta-style" consensus per Gerasimavicius 2023 |
| Evolutionary constraint (MSA-based) | EVcouplings / EVmutation | When Neff ≥ 500 |
| Evolutionary fitness (PLM, MSA-free) | ESM-2 LLR | Ranked #1 of 55 VEPs in Livesey & Marsh 2023 |
| Active-site classification | DSSP + active-site residue lookup | Surface / buried / active-site sphere |

### 3.7 Phenotype × Genotype Correlation (**MISSING FROM CARNIVORENZYME — ADD**)

This is the analytical step that converts a structural atlas into an evolutionary claim.

| Step | Standard tool | Notes |
|------|---------------|-------|
| Phylogenetic generalized least squares (PGLS) | R `caper`, `phylolm` | Controls for phylogenetic non-independence |
| Continuous trait mapping | `phytools` | Visualizes phenotype on tree |
| BayesTraits | independent | Continuous Markov chain analysis |
| Quantitative phenotype | Pitcher pH (Pavlovič 2025) | Mapped across species in species.yaml |

### 3.8 Experimental Validation (Archetype A only)

| Step | Standard tool | Notes |
|------|---------------|-------|
| Gene synthesis | Twist / IDT | $50–$200 per gene |
| Expression | *E. coli* BL21, *Pichia*, baculovirus | Brodelius lab neprosin protocol; Ting 2026 Pichia protocol |
| Kinetic assay | Substrate-specific | Prolyl-endopeptidase activity for neprosin; chitinase activity for GH19 |
| Tm measurement | DSC, DSF | Optional |
| Crystal structure | X-ray, cryo-EM | Not feasible for atlas-scale |

---

## 4. Mapping CarnivorEnzyme Methods to the Standard Stack

| Standard stack element | CarnivorEnzyme current state | Verdict |
|------------------------|------------------------------|---------|
| Sequence retrieval | ✅ fetch_sequences.py (84 verified accessions) | Standard ✓ |
| Genome annotation for un-annotated | ✅ BRAKER2 scripts in workflow/scripts/braker2/ | Standard ✓ |
| MAFFT + trimAl + IQ-TREE | ✅ Phase 2 rules complete | Standard ✓ |
| Ancestral reconstruction | ✅ extract_ancestor.py (IQ-TREE .state) | Standard ✓ |
| **dN/dS selection tests** | ❌ **NOT IMPLEMENTED** | **MISSING — standard requirement** |
| Convergence (site counting) | ✅ detect_convergence.py | Standard ✓ |
| Convergence (ωC error-corrected) | ❌ Cited but not implemented | Optional upgrade |
| Convergence (PLM embedding ACEP) | ❌ Not implemented | Optional |
| AlphaFold3 prediction | ⬜ stub (predict_af3.py) | Standard, needs implementing |
| Ancestral structure prediction | ✅ Implemented (Phase 4A) | Above standard — most papers don't do this |
| FoldX ΔΔG | ⬜ stub | Standard, needs implementing |
| ML ΔΔG (SPURS or RaSP) | ⬜ stub | Above standard — Foldetta-style pair |
| EVcouplings / EVmutation | ⬜ stub | Standard for deep-MSA families |
| ESM-2 / ProSST LLR | ✅ Implemented | Standard ✓ |
| Active-site classification | ⬜ stub (classify_positions.py) | Standard |
| Docking (Vina) | ⬜ stub | Above standard — Butts 2016 precedent only |
| Electrostatics (APBS) | ⬜ stub | Above standard |
| **PGLS phylogenetic comparative** | ❌ **NOT IMPLEMENTED** | **MISSING — Rubisco 2025 precedent** |
| Constant-pH MD | ✅ Implemented | **Over-engineered** — no enzyme evolution paper does this |
| Alchemical FEP | ✅ Implemented | **Over-engineered** — drug discovery tool |
| QM/MM | ✅ Implemented (gated off) | **Over-engineered** — gate off remains correct |
| Wet-lab kinetic validation | ❌ Not planned | Optional — limits journal target if absent |

**Summary:** Standard backbone is intact (sequences → phylogeny → ASR → structure → mutation effects). **Missing:** dN/dS selection analysis and PGLS phenotype correlation — both expected by reviewers. **Over-engineered:** CpHMD + FEP + QM/MM — no published enzyme evolution paper uses these in this context.

---

## 5. The Specific Question CarnivorEnzyme Should Answer

The current framing ("structural atlas of convergent digestive enzymes") is a resource paper. Resource papers exist but are not what the methods stack supports. The right specific question is:

> **"Across nine independent origins of carnivory, do convergent amino acid substitutions in digestive enzymes track adaptation to pitcher fluid pH, and which evolutionary trajectory (functional gain, stability-function tradeoff, or neutral drift) characterizes each substitution?"**

This question:

1. **Is grounded in Fukushima Q1 + Q3a + Pavlovič 2025 Q4.** It directly addresses the open functional question by giving it a quantitative phenotype (pitcher pH).
2. **Is testable with the existing methods stack** (FoldX at multiple pH + PGLS + ancestral reconstruction + convergence detection).
3. **Predicts specific outcomes** (which substitutions matter, which trajectory dominates per family).
4. **Has a wet-lab validation path** (express one ancestral + modern neprosin, measure activity at pH 2.5 vs 7 — single experiment).
5. **Justifies the multi-family scope** — pH adaptation generalizes across families; testing whether the *pattern* repeats across 5 families is a genuinely novel claim.
6. **Maps cleanly to the FoldX × evolutionary quadrant** (the project's core novel analysis).
7. **Reduces over-engineering**: CpHMD becomes optional (supports the pH claim but isn't required); FEP becomes optional (binding affinity is secondary to pH stability).

### Sub-questions the manuscript would answer

| Sub-question | Method | Figure |
|--------------|--------|--------|
| Where are the convergent substitutions located on the structure? | Active-site classification + structural superposition | Fig. 2 |
| Do convergent substitutions destabilize the protein at neutral pH but stabilize it at acid pH? | FoldX ΔΔG at pH 2.5/3.5/5.0 | Fig. 3 |
| Are convergent substitutions evolutionarily favored? | ESM-2 + EVmutation LLR | Fig. 4 |
| Stability vs evolution — what trajectory? | FoldX × EVmutation quadrant | **Fig. 5 (core)** |
| Does the convergent substitution pattern correlate with species' pitcher pH? | PGLS regression | Fig. 6 |
| Are these substitutions under positive selection? | PAML branch-site / HyPhy aBSREL | Fig. 7 |
| Do the convergent substitutions appear at the same point in evolutionary time as carnivory? | dN/dS over time / ancestral resurrection | Fig. 8 |
| Does the ancestral non-carnivorous enzyme function at pH 2.5? *(wet lab option)* | Express ancestor + modern, measure activity | Fig. 9 |

---

## 6. What "Real Outcome" Means

A "real outcome" is a specific, testable claim that the manuscript abstract can make. Three candidate claims, ranked by ambition:

### Outcome A (most defensible, MBE / GBE landing)

> "We show that in **four of five carnivorous plant digestive enzyme families**, convergent amino acid substitutions are biased toward positions that decrease FoldX-predicted ΔΔG at pH 2.5 but not at pH 7, and that this pH-conditional stability shift correlates with species' pitcher fluid pH (PGLS, p < 0.01) — consistent with adaptive evolution for acidic secretion environments."

This is the **default** outcome. Computational-only. Two new analyses required (dN/dS, PGLS). Drop CpHMD/FEP. Achievable by 2026-08-15.

### Outcome B (medium ambition, MBE / NEE)

> Outcome A + "**Resurrection of the predicted MRCA neprosin and a modern *N. gracilis* variant shows the modern enzyme retains >50% activity at pH 2.5 while the ancestor is inactive below pH 4** (kcat/Km, n = 3) — providing the first experimental validation that convergent substitutions in carnivorous plant enzymes drive a quantitative pH-adaptation phenotype."

This adds **one wet-lab experiment** (single ancestor + single modern neprosin). Adds 4–6 months. Journal tier jumps to *MBE* (lock) or *NEE* (probable).

### Outcome C (high ambition, NEE / eLife)

> Outcome B + "**Across the four convergent-substitution families, the proportion of stability-function tradeoff sites (FoldX-destabilizing, evolutionarily-favored) scales with the magnitude of pH adaptation between the ancestor and the most acidic modern species** — suggesting carnivorous plant digestive enzymes evolve through a recurrent stability-cost trajectory rather than through neutral drift or functional gain."

This makes a **trajectory claim** across families. Requires the FoldX × EVmutation quadrant to actually show this pattern (it might not). Higher reward, higher risk. Adds 1–2 wet-lab figures.

---

## 7. The Restructured Paper Outline

Mirroring published Archetype B papers (Fukushima 2017, Rubisco bioRxiv 2025, Wang/Stiffler 2025 PNAS):

### Introduction

- Convergent evolution of carnivory + Fukushima 2017 sequence-level convergence
- Open question (Q3a from Albert/Fukushima 2026)
- Acidic pitcher pH as quantitative phenotype (Pavlovič 2025)
- Specific hypothesis being tested

### Results

1. **Convergent substitutions across 5 families × 9 origins** (Fig. 1: phylogeny with convergent sites colored)
2. **Structural mapping** (Fig. 2: AF3 structures with convergent positions; active-site/surface/buried classification)
3. **pH-conditional stability shifts** (Fig. 3: FoldX ΔΔG at pH 2.5/3.5/5.0)
4. **Evolutionary constraint** (Fig. 4: EVmutation + ESM-2 LLR distributions)
5. **Stability-evolution trajectory** (Fig. 5: FoldX × EVmutation quadrant — primary novel figure)
6. **Phylogenetic correlation with pitcher pH** (Fig. 6: PGLS regression)
7. **Positive selection signatures** (Fig. 7: PAML branch-site, HyPhy aBSREL)
8. **Optional: ancestral resurrection** (Fig. 8: kcat/Km of MRCA vs modern at pH 2.5)

### Discussion

- Comparison to Fukushima 2017, Rubisco convergence, Wang/Stiffler 2025
- Mechanism: stability-function tradeoff as the dominant trajectory
- Limitations (computational predictions, one wet-lab validation, MSA depth caveats)
- Future work: CpHMD/FEP/QM-MM as preliminary data for mechanism follow-up

### Methods

- Sequence retrieval (NCBI / ENA / JGI)
- Phylogeny (MAFFT + trimAl + IQ-TREE 3)
- Ancestral reconstruction (IQ-TREE marginal)
- Convergence detection (Fukushima 2017 method + ωC + ACEP)
- Selection tests (PAML, HyPhy)
- AF3 structure prediction
- FoldX 5.1 ΔΔG + SPURS or RaSP (Foldetta consensus)
- EVmutation + ESM-2 LLR
- PGLS comparative analysis
- (Wet lab if included)

### Supplementary

- CpHMD pKa shifts (mechanism preliminary data)
- FEP ΔΔG_bind (substrate affinity preliminary data)
- QM/MM transition state (top-1 sites only)
- Docking + electrostatics
- Expression quantification

---

## 8. Sources

- Fukushima K et al. (2017) *Nat. Ecol. Evol.* 1:0059
- Fukushima K, Pollock DD. (2023) *Nat. Ecol. Evol.* 7:1156
- Albert VA, Avila-Robledillo L, Fleck S, et al. (2026) *Trends in Genetics* doi:S0168-9525(26)00034-X
- Pavlovič A (2025) *New Phytologist* doi:10.1111/nph.70229
- Freund M et al. (2022) *Plant Physiology* 190:44
- Iñiguez C et al. (2022) *Sci. Adv.* doi:10.1126/sciadv.abm6871
- Wang Z et al. (2025) *PNAS* doi:10.1073/pnas.2418254122 (ACEP)
- Rubisco convergence (2025) bioRxiv 2025.10.08.681247
- Butts CT et al. (2016) *Proteins* 84:1517
- Riziotis IG et al. (2025) *FEBS Journal* doi:10.1111/febs.17332
- Gerasimavicius L, Livesey BJ, Marsh JA (2023) *Protein Science* doi:10.1002/pro.4688 (Foldetta)
- McKeown AN, Thornton JW (2024) *Nature* (DMS on resurrected ancestors)
- JACS Au (2024) doi:10.1021/jacsau.4c00653 (ASR review)
- PAML codeml (Yang 2007)
- HyPhy aBSREL (Smith et al. 2015)
- phytools R package (Revell 2024)

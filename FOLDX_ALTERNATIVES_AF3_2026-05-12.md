# FoldX Alternatives for AF3-Predicted Structures — Targeted Literature Review

> Date: 2026-05-12
> Follows: [FOLDX_REVIEW_2026-05-12.md](FOLDX_REVIEW_2026-05-12.md)
> Question asked: Is there a ΔΔG predictor that works *better than FoldX specifically on AF3-predicted structures* AND handles pH dependency?

---

## 1. Direct Answer

**No single method published as of May 2026 simultaneously (a) is trained or validated on AlphaFold-predicted inputs, (b) handles pH-dependent ΔΔG, AND (c) outperforms FoldX 5.1 on heterogeneous plant enzymes.** Each candidate gives up at least one of those three.

But the search surfaced **three specific methods** that should be added to CarnivorEnzyme's stability stack, plus **one critical reality check** about AF3 native confidence as a ΔΔG proxy.

---

## 2. The Three New Methods to Add

### 2.1 SPURS — the AF-native ΔΔG predictor
**Citation:** Cao et al. (2025) Generalizable and scalable protein stability prediction with SPURS. *Nature Communications* 16. doi:10.1038/s41467-025-67609-4
**Preprint:** bioRxiv 2025.02.13.638154
**Code:** github.com/luo-group/SPURS

**Why it matters:** **SPURS is the only published ΔΔG method explicitly designed and trained on AlphaFold-predicted backbones.** Rewires ProteinMPNN structural priors into ESM-2 embeddings. The Cao et al. paper validates on AF-predicted inputs as the primary use case, not as an afterthought.

**Performance:**
- SCC = 0.54 on Human Domainome (best baseline ThermoMPNN: 0.49)
- Outperforms ThermoMPNN / RaSP / Stability Oracle on 7 of 8 standard test sets
- Open source, pip-installable

**What it lacks:** No pH input. Trained on cDNA-display proteolysis at near-neutral pH (Tsuboyama 2023 Megascale conditions).

**Verdict:** Best single deep-learning alternative if AF3 input quality is the limiting concern. Replaces or complements ThermoMPNN in the planned ensemble.

### 2.2 CatOpt — closest published precedent to the biological question
**Citation:** Wang et al. (2025) CatOpt: Sequence-based design of acid-resistant enzyme mutations. *ACS Synthetic Biology*. doi:10.1021/acssynbio.5c00679

**Why it matters:** Directly addresses CarnivorEnzyme's biological question — given an enzyme sequence, predict pHopt and design single-point mutations that shift it acidward. **Experimentally validated** by engineering a 7% activity improvement at low pH in a diacetylchitobiose deacetylase.

**Performance:**
- R² = 0.479 on pHopt prediction hold-out
- Validated wet-lab confirmation in chitin-related enzyme — very close functional analogy to GH19 chitinases

**What it outputs:** Not ΔΔG. Predicts shift in pHopt upon mutation — arguably more biologically meaningful than ΔΔG_folding for digestive enzymes that must function at pH 2.5.

**Verdict:** Add as **fourth axis** alongside FoldX × SPURS × ProSST. The "does this mutation shift pHopt acidward?" signal is what convergent substitutions in carnivorous plant enzymes likely encode, and no other tool measures this directly.

### 2.3 PypKa — better upstream protonation than PROPKA
**Citation:** Reis & Machuqueiro (2024) PypKa server: pKa predictions and structural protonation across AlphaFoldDB. *Nucleic Acids Research* 52:W294. doi:10.1093/nar/gkae255

**Why it matters:** PROPKA (currently planned for upstream of FoldX) has ~1 pH-unit prediction error and is structure-quality sensitive. PypKa was **explicitly benchmarked on AlphaFold structures**, has precomputed pKa values across the AlphaFoldDB, and offers ML acceleration. Combines Poisson-Boltzmann electrostatics with ML.

**Performance:**
- Top non-MD method in Reis 2024 benchmark
- Has precomputed pKa for >200M proteins in AlphaFoldDB
- ML mode is ~100× faster than full PB calculation

**Verdict:** **Swap PROPKA → PypKa** upstream of FoldX scan at each pH. Field-aware best practice for AF inputs.

---

## 3. The Critical Reality Check: AF3 Confidence is NOT a Stability Proxy

A bioRxiv preprint (May 2024.05.25.595871, "AlphaFold3, a secret sauce for predicting mutational effects on protein-protein interactions") and follow-up **Zhang et al. *JCIM* 2024** (doi:10.1021/acs.jcim.4c00976) report PCC = 0.86 for AF3-derived ΔΔG_binding on **protein-protein interfaces**.

**This does not generalize to monomer folding stability.**

- Akdel et al. 2022 (NSMB) explicitly showed that **ΔpLDDT does NOT correlate with experimental ΔΔG** for monomer folding mutations.
- Pak et al. 2023 (*PLOS ONE*) confirmed the same finding.
- No published paper claims AF3 internally predicts monomer ΔΔG of folding.

**Implication:** Do not attempt to use AF3 pLDDT or PAE shifts as a stability signal for convergent substitution analysis. Use AF3 confidence **only as a filter** (exclude pLDDT < 70 sites from quantitative ΔΔG claims), not as a predictor.

---

## 4. Complete Tier Ranking (Refined)

### Tier A — AF-compatible structure-based ΔΔG (neutral pH)

| Rank | Method | AF-validated | PCC range | Why |
|------|--------|--------------|-----------|-----|
| 1 | **SPURS** (Cao 2025, *Nat. Commun.*) | **explicit training on AF** | 0.54 (Domainome) | AF-native by design |
| 2 | RaSP (Blaabjerg 2023, *eLife*) | yes (0.94 vs 0.97 X-ray) | 0.71 (Megascale) | Fastest; saturation maps |
| 3 | Stability Oracle (Diaz 2024, *Nat. Commun.*) | yes (ColabFold tested) | ~0.55 | Graph transformer |
| 4 | ThermoMPNN-D (Dieckhaus 2024, *PNAS*) | transfer (Akdel-style) | 0.54–0.59 | Strong baseline |
| 5 | Cagiada 2025 (*Protein Science*) | AF/ESMFold native | R=0.7 (absolute) | Predicts absolute ΔG |

### Tier B — Sequence-only (structure quality irrelevant)

| Method | Citation | Use |
|--------|----------|-----|
| JanusDDG | Ressa et al. 2026, *Commun. Biol.* (arXiv 2504.03278) | Antisymmetric, transitive; works on low-identity test sets |
| PROSTATA | Umerenkov 2023 | Transformer on sequence |
| DDMut / DDMut-PPI | Zhou 2023/2024, *NAR* | Web only; explicit pH 7 limitation |

### Tier C — Plant / acid-pH specific (the new axis)

| Method | Citation | Output | Acid-pH relevance |
|--------|----------|--------|-------------------|
| **CatOpt** | Wang 2025, *ACS Synth. Biol.* | Δ pHopt per mutation | **Direct — wet-lab validated** |
| EpHod | Gado 2023, bioRxiv | pHopt prediction | Validates ortholog pHopt assignments |
| TemStaPro | Pudžiuvelytė 2024, *Bioinformatics* | Tm class | Global Tm sanity check |
| DeepSTABp | Jung 2023, *NAR* | Tm regression | Global Tm sanity check |
| DeepTM | Liu 2023, *CSBJ* | Tm of plant/PET enzymes | Best on PET-degrading enzymes |

### Tier D — Explicit pH-dependent ΔΔG

| Method | Citation | Notes |
|--------|----------|-------|
| **FoldX 5.1** | Botte 2025, *Bioinformatics* 41:btaf064 | **The only one. RETAIN.** |
| Nonequilibrium alchemy pKa | Wilson 2023, *JCTC* | Most accurate; ~1000× more expensive |
| GROMACS CpHMD | Gapsys 2022, *JCTC* | Already in Phase 5E |

### Tier E — Native AF / DeepMind heads

| Method | Citation | Verdict |
|--------|----------|---------|
| AlphaMissense | Cheng 2023, *Science* | DeepMind explicitly says it doesn't predict stability |
| AF3 ΔPAE for monomer ΔΔG | various | **Does not work** per Akdel 2022 |
| AF3 ΔPAE for PPI ΔΔG_bind | Zhang 2024, *JCIM* | PCC 0.86 for **binding** — not folding |
| AlphaProteo | DeepMind 2024, *Nature* | Binder design only |

---

## 5. The Recommended Final Stack

```
Stability axes (4-method consensus):
  1. FoldX 5.1 at pH 2.5 / 3.5 / 5.0     ← pH-aware physics (KEEP)
       └─ upstream: PypKa for protonation (UPGRADE from PROPKA)
  2. SPURS                                 ← AF-native deep learning (ADD)
  3. RaSP                                  ← saturation maps, cheap (ADD)
  4. CatOpt                                ← pHopt shift, biological signal (ADD)
       sanity check: DDGun3D               ← untrained baseline

Evolutionary axes (3-method consensus):
  1. ProSST                                ← structure-aware PLM (primary)
  2. SaProt                                ← Foldseek 3Di (secondary)
  3. EVE / EVcouplings                     ← deep-MSA families only

Active-site protonation:
  CpHMD (Phase 5E)                         ← Glu/Asp pKa at acid pH
  Nonequilibrium alchemy pKa                ← top-5 FEP targets only
```

### Why this stack is reviewer-defensible

1. **FoldX retains the pH moat.** No DL method handles pH; reviewers cannot complain about FoldX as the pH axis if no alternative exists.
2. **SPURS is explicitly AF-trained.** Pre-empts "why didn't you use a 2024-era DL method on AF inputs?" — SPURS is the 2024-era DL method designed for AF inputs.
3. **CatOpt is wet-lab-validated on a closely-related chitinase substrate.** Directly addresses pH adaptation, the biological question.
4. **The 4-axis consensus is novel.** No published carnivorous plant or general plant enzyme convergence study uses 4 stability axes.

---

## 6. Concrete Implementation Tasks

| Task | File | Effort |
|------|------|--------|
| Tighten `plddt_exclude: 50 → 70` | [config/config.yaml](config/config.yaml) | 1 line |
| Add `score_thermompnn.py` (now Tier 2) | [workflow/scripts/](workflow/scripts/) | ~1 day |
| Add **`score_spurs.py`** (PRIMARY DL ΔΔG) | [workflow/scripts/](workflow/scripts/) | ~1 day |
| Add `score_rasp.py` | [workflow/scripts/](workflow/scripts/) | ~1 day |
| Add **`score_catopt.py`** (Δ pHopt) | [workflow/scripts/](workflow/scripts/) | ~half day |
| Swap PROPKA → **PypKa** upstream of FoldX | [run_foldx_scan.py](workflow/scripts/run_foldx_scan.py) (when implementing) | ~2 hours |
| New env file `workflow/envs/stability_ensemble.yaml` (SPURS + RaSP + CatOpt + DDGun3D) | new file | ~1 hour |
| New Snakemake rule `stability_ensemble.smk` | [workflow/rules/](workflow/rules/) | ~half day |
| Rename `compare_foldx_prosst.py` → `compare_stability_evolution.py` with 4 stability + 3 evolution axes | renamed | ~1 day |
| Update [TASKS.md](TASKS.md) Tier 1 with these scripts | TASKS.md | ~5 min |

---

## 7. What This Changes for the Manuscript

The original CarnivorEnzyme novelty pitch (per [LITERATURE_REVIEW_2026-05-12.md](LITERATURE_REVIEW_2026-05-12.md)) was:

> "FoldX × evolutionary-constraint axes on convergent substitutions across multi-family multi-lineage data."

The revised pitch becomes:

> "**Four-axis stability ensemble × three-axis evolutionary ensemble** on convergent substitutions, integrating pH-aware physics (FoldX 5.1), AF-native deep learning (SPURS), saturation-scale stability landscape (RaSP), and biologically-direct acid-adaptation prediction (CatOpt), against structure-aware (ProSST), Foldseek-aware (SaProt), and DCA-based (EVE/EVcouplings) evolutionary constraint."

The expanded methods stack is more defensible to reviewers, more biologically resolved, and more genuinely novel than any single-axis or two-axis approach in the literature.

---

## 8. Open Questions for Next Search

1. Does AF3 export per-residue confidence (PAE/pLDDT) compatible with SPURS input expectations? Need to verify SPURS accepts ColabFold/AF3 outputs natively, not just AF2.
2. Is CatOpt API available, or does it require local install of a heavy model? Need to check `https://github.com/microsoft/CatOpt` or supplementary code.
3. Has any 2026 preprint extended SPURS to pH-conditional training? (Worth setting up a Google Scholar alert.)

---

## 9. Sources Added Beyond FOLDX_REVIEW_2026-05-12.md

- Cao et al. (2025) SPURS: Generalizable and scalable protein stability prediction. *Nat. Commun.* 16. doi:10.1038/s41467-025-67609-4 (bioRxiv 2025.02.13.638154)
- Wang et al. (2025) CatOpt: acid-resistance enzyme design. *ACS Synth. Biol.* doi:10.1021/acssynbio.5c00679
- Reis & Machuqueiro (2024) PypKa server. *NAR* 52:W294. doi:10.1093/nar/gkae255
- Ressa et al. (2026) JanusDDG thermodynamic consistency. *Commun. Biol.* doi:10.1038/s42003-026-09632-9 (arXiv 2504.03278)
- Zhang et al. (2024) AlphaFold3 PPI ΔΔG. *JCIM*. doi:10.1021/acs.jcim.4c00976
- Gado et al. (2023) EpHod pHopt prediction. bioRxiv 2023.06.22.544776
- Liu et al. (2023) DeepTM thermostability of PET-degrading enzymes. *CSBJ*
- Pudžiuvelytė et al. (2024) TemStaPro. *Bioinformatics* 40:btae157
- Jung et al. (2023) DeepSTABp. *NAR*

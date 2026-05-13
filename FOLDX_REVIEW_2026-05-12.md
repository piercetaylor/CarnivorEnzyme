# FoldX vs Alternatives — Literature Findings

> Date: 2026-05-12
> Question asked: Is FoldX 5.1 still state-of-the-art? Is FoldX-on-AF3 standard practice? Are there better pH-aware alternatives?
> Source: Strategic literature search May 2026 (see [LITERATURE_REVIEW_2026-05-12.md](LITERATURE_REVIEW_2026-05-12.md))

---

## 1. Top-Line Findings

| Question | Verdict |
|----------|---------|
| Is FoldX 5.1 still SOTA for ΔΔG? | **No.** Deep-learning methods (ThermoMPNN, SPURS, PROSTATA) outperform FoldX by PCC 0.30+ on Megascale and S669 benchmarks. |
| Is FoldX-on-AF3 publishable? | **Yes**, with explicit caveats — Akdel 2022 (*NSMB*) is the standard citation. |
| Is FoldX's pH dependency useful for pitcher fluid pH? | **Yes** — it's the only major method with pH awareness, but Botte 2025 validation is thin at acid pH. |
| Should CarnivorEnzyme replace FoldX? | **No** — add deep-learning axes alongside it (ensemble approach). |

---

## 2. Benchmark Numbers (Megascale + S669 + FireProt HF)

The Tsuboyama et al. 2023 Megascale dataset (*Nature* 620:434, ~776,000 ΔΔG measurements across 479 protein domains, cDNA Display Proteolysis) reshuffled the field. Independent benchmarks as of May 2026:

| Method | Megascale PCC | S669 PCC | Speed | Type |
|--------|---------------|----------|-------|------|
| **ThermoMPNN** (Dieckhaus 2024, *PNAS*) | ~0.75 | ~0.43 | fast | ProteinMPNN backbone + fine-tune |
| **SPURS** (Cai 2025, *Nat. Commun.*) | beats above on 7/8 sets | best on S669 | fast | ESM2 + ProteinMPNN rewiring |
| **PROSTATA** (Umerenkov 2023, *Bioinformatics*) | ~0.72 | ~0.51 | fast | Transformer |
| **RaSP** (Blaabjerg 2023, *eLife*) | ~0.71 | ~0.39 | 480–1000× faster than Rosetta | DL |
| **Pythia** (2024, *The Innovation*) | — | Spearman 0.66 | ultrafast | self-supervised |
| **DDMut** (Zhou 2023, *NAR*) | — | ~0.70 (own test) | web only | siamese DL |
| Rosetta cartesian_ddg (Park 2016) | ~0.53 | ~0.50 | slow | physics |
| **FoldX 5.0** | ~0.40 | ~0.30 (S461) | slow (~1 mut/min) | physics |
| ThermoNet (Li 2020, *PNAS*) | ~0.33 | — | medium | 3D CNN |
| DDGun3D (Montanucci 2019) | — | ~0.43 | fast | untrained, robust |
| Stability Oracle (Diaz 2024, *Nat. Commun.*) | comparable to ThermoMPNN | — | medium | transformer |
| Cagiada 2025 (absolute ΔG, *Protein Science*) | — | — | — | predicts absolute ΔG, MAE 1.5 kcal/mol |

**Botte 2025 revision** improved FoldX 0.693 → 0.711 PCC on its own held-out FSD test (real but modest gain). Does **not** close the gap with ThermoMPNN/SPURS on Megascale or S669.

**Buß 2018 audit numbers** (29% success for stabilizing, 69% for destabilizing) have not been re-validated for FoldX 5.1.

---

## 3. FoldX-on-AlphaFold Is Now Standard Practice

The key published assessment: **Akdel et al. 2022** (*Nat. Struct. Mol. Biol.* 29:1056, doi:10.1038/s41594-022-00849-w). Community assessment across >100,000 DMS mutations and three predictors (FoldX, Rosetta, DynaMut2). **AF2-derived ΔΔG matched or exceeded experimental-structure predictions.**

Supporting work:

- **Marabotti et al. 2022** (*Bioinformatics* 38:4312) — AF models outperform homology models for FoldX in most cases.
- **Frontiers in Genetics 2023** (doi:10.3389/fgene.2023.1052383) — explicitly benchmarked FoldX-on-AF for cancer missense variants; usable performance.

### Documented pitfalls reviewers will probe

1. **pLDDT filtering is essential.** Field-standard is **pLDDT ≥ 70** at the mutated residue. CarnivorEnzyme [config/config.yaml](config/config.yaml) currently has `plddt_exclude: 50` — too permissive for quantitative ΔΔG claims. **ACTION: tighten to 70.**
2. **RepairPDB on AF inputs.** Five rounds (current setting) is correct; some labs run up to 10 for AF.
3. **Pak et al. 2023** (*PLOS ONE*, doi:10.1371/journal.pone.0282689) — warned against using AF's own pLDDT/ΔpLDDT shifts as a stability proxy. **Does not apply** to FoldX-on-AF since we're running real FoldX, not pLDDT differential.
4. **Buel & Walters 2022** (*Nat. Struct. Mol. Biol.* 29:1) — AF structures are conformational averages. Mitigation: 5-seed AF3 protocol generates a structural ensemble.
5. **Tsuboyama 2023 used X-ray structures**, not AF. All Megascale benchmarks are upper bounds for AF input.

---

## 4. The pH Dependency — Where FoldX Still Uniquely Wins

**No deep-learning ΔΔG method handles pH-dependent protonation as of May 2026.** ThermoMPNN, SPURS, RaSP, PROSTATA, ProteinMPNN are all trained on neutral-pH data with no input for pH state.

For carnivorous plant enzymes at pitcher fluid pH 1.9–5.0, this is decisive:

- The whole project question is intrinsically pH-dependent.
- FoldX 5.1 (Botte 2025) added experimentally-measured pKas for Asp/Glu/Lys/Arg/Tyr (previously only His/Cys).
- **Botte 2025 validates on FSD at pH 6.5–8.0 — not benchmarked at pH 2.5–5.0.** So acid-pH FoldX numbers carry more uncertainty than neutral.

### Upstream protonation tools (use before FoldX)

| Tool | Citation | Use |
|------|----------|-----|
| **PypKa** | Reis et al. 2024, *NAR* 52:W294 | Poisson-Boltzmann pKa across AlphaFoldDB; top non-MD method |
| **PROPKA 3.4** (via PDB2PQR) | Olsson 2011, *JCTC* 7:525 | Long-standing standard; ~1 pH unit error |
| **DeepKa, pKa-ANI** | Reis 2024 benchmark | ML pKa predictors; pKa-ANI top performer |
| **Constant-pH MD** (GROMACS λ-dynamics) | Gapsys et al. 2022, *JCTC* 18:6320 | Gold standard; computationally expensive |
| **Nonequilibrium alchemy pKa** | Wilson 2023, *JCTC* doi:10.1021/acs.jctc.3c00721 | Most accurate, expensive |

### CpHMD validation for aspartic proteases

**Aguilar et al. 2010** (PMC7312390) — explicitly validated CpHMD for aspartic protease catalytic-site pKa. Directly relevant to nepenthesins/droserasins (Asp32/Asp215 numbering) and neprosin (E188/E297). **Phase 5E [run_cphmd.py](workflow/scripts/run_cphmd.py) implementation already correct.**

---

## 5. The Recommended Ensemble Approach

**Don't replace FoldX. Add deep-learning axes alongside it.**

### Stability axes (new ensemble)

1. **FoldX 5.1** at pH 2.5 / 3.5 / 5.0 — pH-aware physics (RETAIN)
2. **ThermoMPNN** (or SPURS) — DL state-of-art (ADD)
3. **RaSP** — saturation scanning, cheap (ADD)
4. **DDGun3D** — untrained baseline, robust to AF quality (ADD as sanity check)

### Evolutionary axes (per [LITERATURE_REVIEW_2026-05-12.md](LITERATURE_REVIEW_2026-05-12.md))

1. **ProSST** — structure-aware PLM (primary)
2. **SaProt** — Foldseek 3Di-aware (secondary)
3. **EVE** / **EVcouplings** — deep-MSA families only

### Active-site protonation

- **CpHMD** (Phase 5E, existing) — Glu/Asp pKa at acid pH

### What this enables

The original FoldX × ProSST 2D quadrant becomes a **multi-axis classification** of each convergent site:

- Consensus stability cost (FoldX + ThermoMPNN + RaSP)
- pH sensitivity (FoldX ΔΔG pH 2.5 vs pH 5.0; CpHMD pKa shift)
- Evolutionary constraint (ProSST + SaProt + EVE LLR)
- Substrate affinity (FEP, top ≤20 sites)

A site destabilizing on all stability axes, evolutionarily favored on all PLMs, and pKa-shifted toward pitcher pH is a **very high-confidence adaptive convergent substitution**. Stronger than any single-method claim.

---

## 6. Concrete Pipeline Changes

| Change | File | Effort |
|--------|------|--------|
| Tighten `plddt_exclude: 50 → 70` | [config/config.yaml](config/config.yaml) | 1 line |
| Add `score_thermompnn.py` | [workflow/scripts/](workflow/scripts/) | ~1 day |
| Add `score_rasp.py` | [workflow/scripts/](workflow/scripts/) | ~1 day |
| Add `score_ddgun3d.py` (baseline) | [workflow/scripts/](workflow/scripts/) | ~half day |
| Add upstream PypKa/PROPKA call before FoldX scan | `run_foldx_scan.py` (when implementing) | ~2 hours |
| Extend `compare_foldx_prosst.py` → `compare_stability_evolution.py` with 4 stability axes | renamed | ~1 day |
| Update [TASKS.md](TASKS.md) with new scripts in Tier 1 | TASKS.md | ~5 min |
| Cite Akdel 2022, Pak 2023, Buel & Walters 2022, Tsuboyama 2023 in Methods | manuscript draft | future |

---

## 7. Reviewer Defensibility

This ensemble approach has not been published for any plant enzyme convergence study. It pre-empts the predictable reviewer challenge: *"Why didn't you use a 2024-era deep-learning method?"*

Defensible for *Nature Communications*, *Molecular Biology and Evolution*, *Nature Ecology & Evolution* if:

1. Akdel 2022 cited as FoldX-on-AF justification
2. pLDDT ≥ 70 enforced at convergent sites
3. 7ZVA crystal-vs-AF3 benchmark passes r > 0.65 (Test Gate 4 already does this)
4. Report ΔΔG distributions and rankings — not absolute kcal/mol claims
5. At least one DL stability predictor present as a second axis
6. FoldX retained as the pH-aware physics axis

---

## 8. Follow-Up Search (Completed 2026-05-12)

**Question:** Is there a single ΔΔG predictor that beats FoldX *specifically on AF3-predicted structures* AND handles pH dependency?

**Answer:** No single tool wins on all three criteria. But three new methods should be added to the stack — see [FOLDX_ALTERNATIVES_AF3_2026-05-12.md](FOLDX_ALTERNATIVES_AF3_2026-05-12.md) for full analysis.

Top finding: **SPURS** (Cao et al. 2025, *Nat. Commun.*) is the only published ΔΔG method explicitly trained on AlphaFold-predicted backbones. SCC 0.54 on Human Domainome, beats ThermoMPNN (0.49). Should be added as the primary deep-learning axis.

Second finding: **CatOpt** (Wang et al. 2025, *ACS Synth. Biol.*) is wet-lab-validated for designing acid-resistant mutations in a chitinase-related enzyme. Directly addresses the biological question of pH adaptation. Add as a 4th stability axis.

Third finding: **PypKa** (Reis & Machuqueiro 2024, *NAR* 52:W294) is the AF-validated upstream protonation tool. Replace planned PROPKA with PypKa.

Critical reality check: **AF3 pLDDT/PAE shifts do NOT correlate with monomer ΔΔG** per Akdel 2022 + Pak 2023. The Zhang 2024 JCIM paper showing PCC 0.86 for AF3-derived ΔΔG_binding applies only to protein-protein interfaces, not monomer folding. Do not use AF3 confidence as a stability proxy.

---

## 9. Sources

- Tsuboyama K et al. (2023) Mega-scale experimental analysis of protein folding stability. *Nature* 620:434. doi:10.1038/s41586-023-06328-6
- Akdel M et al. (2022) A structural biology community assessment of AlphaFold2. *Nat. Struct. Mol. Biol.* 29:1056. doi:10.1038/s41594-022-00849-w
- Botte M et al. (2025) FoldX force field revision with pH dependency. *Bioinformatics* 41:btaf064. doi:10.1093/bioinformatics/btaf064
- Dieckhaus H et al. (2024) ThermoMPNN: Transfer learning for protein stability prediction. *PNAS* doi:10.1073/pnas.2314853121
- Cai et al. (2025) SPURS: Stability prediction using rewiring. *Nat. Commun.* doi:10.1038/s41467-025-67609-4
- Blaabjerg LM et al. (2023) Rapid protein stability prediction using deep learning (RaSP). *eLife* doi:10.7554/eLife.82593
- Diaz DJ et al. (2024) Stability Oracle. *Nat. Commun.* doi:10.1038/s41467-024-49780-2
- Umerenkov D et al. (2023) PROSTATA. *Bioinformatics* doi:10.1093/bioinformatics/btad671
- Zhou Y et al. (2023) DDMut. *NAR* doi:10.1093/nar/gkad472
- Pancotti C et al. (2022) Best-templates vs homology models for FoldX. *Bioinformatics* 38:4312
- Pak MA et al. (2023) Using AlphaFold for stability assessment caveats. *PLOS ONE* doi:10.1371/journal.pone.0282689
- Buel GR, Walters KJ (2022) Limitations of AlphaFold. *Nat. Struct. Mol. Biol.* 29:1
- Cagiada M et al. (2025) Absolute ΔG prediction. *Protein Science* doi:10.1002/pro.5233
- Reis P et al. (2024) PypKa server. *NAR* 52:W294. doi:10.1093/nar/gkae412
- Aguilar B et al. (2010) Catalytic-site pKa via CpHMD on aspartic proteases. PMC7312390
- Buß O et al. (2018) FoldX as protein engineering tool. *Comput. Struct. Biotechnol. J.* 16:25

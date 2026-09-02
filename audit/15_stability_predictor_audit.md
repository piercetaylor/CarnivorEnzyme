# 15 — Stability predictor audit: what to keep, what to cut, what was never sound

> Date: 2026-09-01
> Scope: the ΔΔG / stability axis of the CarnivorEnzyme method stack — FoldX, SPURS, RaSP,
> ThermoMPNN, CatOpt, DDGun3D and the rest — plus the pH-dependent stability claim that axis was
> supposed to carry.
> Supersedes: `archive/FOLDX_REVIEW_2026-05-12.md`, `archive/FOLDX_ALTERNATIVES_AF3_2026-05-12.md`,
> and `PROJECT_PLAN.md` §2.2 (the 2026-08-25 "SPURS only, FoldX dropped entirely" decision).
> Method: two independent literature-verification passes resolving every citation against Crossref,
> PubMed/PMC, Europe PMC, OpenAlex, bioRxiv/openRxiv and the GitHub API; plus first-hand FoldX 5.1
> runs on the neprosin crystal structure 7ZVA.

---

## 0. Verdict in one page

| Tool | Prior status in repo | Verdict | Why |
|------|---------------------|---------|-----|
| **FoldX 5.1** | Dropped entirely (2026-08-25) | **REINSTATE as primary — re-scoped** | The only predictor published as retaining correlation at *surface* residues, which is where this project's substitutions are. Already licensed and installed. But: single reference pH, ranked output, no kcal/mol thresholds. |
| **RaSP** | Proposed then dropped | **ADD as the one cross-check** | The only tool with published crystal-vs-AlphaFold concordance (ρ̄ = 0.94 vs. 0.97 crystal-to-crystal). Directly answers the one thing FoldX is weakest at. |
| **SPURS** | Sole stability axis | **DEMOTE to optional third; not primary, not sole** | Real and good, but has no experimental-vs-AF ablation, no pH input, is trained on ≤72-residue domains, and has 15 citations. Sole-axis status is not supportable. |
| **ThermoMPNN, PROSTATA, Pythia, DDMut, JanusDDG, DDGun3D, Stability Oracle** | Proposed as ensemble axes | **CUT** | Adding them buys correlated error, not information (§4.3). |
| **CatOpt** | Proposed as a 4th "acid-adaptation" axis | **CUT — methodologically invalid here** | Whole-sequence pH-optimum regressor with RMSE 0.833 pH units and no per-mutation head. A single substitution cannot move it beyond noise. |
| **EpHod** | Proposed supporting axis | **CUT for per-site use** | Same defect: whole-sequence, RMSE 1.25 pH units. Usable only to sanity-check whole-ortholog pH optima. |
| **PypKa (replacing PROPKA)** | Proposed upgrade | **CUT the swap; keep PROPKA** | PypKa beats PROPKA by 0.04 pH units RMSE and beats assuming solution pKa by 0.02. There is no upgrade here. |
| **FoldX pH sweep (2.5/3.5/5.0)** | Planned headline result (Fig. 3, H1) | **CUT as a primary result** | Measured directly: the pH 2.5 → 5.0 shift at exposed sites averages **0.12 kcal/mol**, ~10× below FoldX's own RMSE, 0/43 mutations clearing it (§1.1). The module's own validation is one histidine in one protein and its increment was not statistically significant (§3.2). |
| **CpHMD (Phase 5E)** | "Sole pH-resolved stability signal" | **DOWNGRADE to exploratory; cannot carry the pH claim as configured** | Runs only on an unreleased beta fork whose authors ask users not to publish from it; 3 pH points cannot produce a pKa; budget underscoped 8–17× (§3.3). |
| **Electrostatics (Phase 6)** | Stub, deprioritised | **PROMOTE to the primary pH axis** | It is what Fukushima 2017's own interpretation and the acidophile-adaptation literature both point to (§5). |

**The single most important finding is not about any tool.** It is that the paper this project is
built on already reported that the convergent substitutions are surface-exposed, are not near
catalytic residues, and act on solution interactions *rather than protein conformation* — which is
the observable ΔΔG-of-folding measures. See §2.

---

## 1. What was verified first-hand

FoldX 5.1 is installed locally (`C:\FoldX\foldx.exe`) and its `--pH`, `--temperature` and
`--ionStrength` parameters are real. To test whether the stability axis carries usable signal at
the site class this project cares about, RepairPDB was run on chain A of **7ZVA** (neprosin zymogen,
1.80 Å), relative solvent accessibility was computed per residue (Shrake–Rupley, normalised by the
Tien et al. 2013 Gly-X-Gly maxima — the same normalisation Fukushima et al. 2017 used), and 16
positions were selected in four classes: buried/exposed × titratable/non-titratable. Each was
mutated to a fixed six-residue panel (A, D, E, K, R, N) and scanned with `--command=PositionScan`.

**FoldX ΔΔG at pH 2.5, 7ZVA chain A, 87 mutations:**

| class | n | mean rel. SASA | range (kcal/mol) | mean \|ΔΔG\| | frac \|ΔΔG\| > 1.0 | frac \|ΔΔG\| > 1.25 | frac \|ΔΔG\| > 2.9 |
|---|---|---|---|---|---|---|---|
| buried, non-titratable | 23 | 0.00 | 3.21 … 14.60 | 7.07 | 1.00 | 1.00 | 1.00 |
| buried, titratable | 21 | 0.00 | 0.82 … 14.49 | 6.10 | 0.95 | 0.95 | 0.91 |
| exposed, non-titratable | 22 | 0.89 | −1.59 … 0.77 | 0.42 | 0.09 | 0.09 | 0.00 |
| exposed, titratable | 21 | 0.69 | −0.72 … 0.72 | 0.39 | 0.00 | 0.00 | 0.00 |

Reference points for the last three columns: **1.0 kcal/mol** is this repo's own
`config.yaml: foldx.ddg_destabilizing` threshold; **1.25 kcal/mol** is FoldX 5.1's own reported
RMSE (Delgado et al. 2025); **2.9 kcal/mol** is the published 95% prediction interval for FoldX
folding ΔΔG (Sapozhnikov et al. 2023).

**All 43 exposed-site mutations lie within a ±1.6 kcal/mol band, mean |ΔΔG| = 0.40 kcal/mol.
Zero of 43 exceed the 95% prediction interval. Two of 43 exceed even the point RMSE.**

The consequence for the pipeline as designed is concrete: against the configured thresholds, every
buried position is classified "destabilizing" and essentially every surface position is classified
"neutral". The stability axis of the quadrant plot is a near-constant at exactly the positions the
project is about. This is not a tuning problem — the tool's entire dynamic range at the surface is
roughly a third of its own error bar.

*(Reproduce: `audit/scripts/` is not yet populated; the scan harness used here is retained in the
session scratchpad. Re-running it against the AF3 predictions once Phase 3 lands is TODO Tier 1.)*

### 1.1 pH sensitivity, measured

The same 87 mutations were then rescanned at **pH 5.0** and **pH 7.0** — the identical structure,
identical positions, only the `--pH` argument changed. This is a direct test of the planned Figure 3
/ H1 result ("do convergent substitutions destabilize at pH 5.0 but not at pH 2.5?").

**ΔΔG(pH 2.5) − ΔΔG(pH 5.0), by class:**

| class | n | mean Δ | mean \|Δ\| | max \|Δ\| | frac \|Δ\| > 1.25 (FoldX RMSE) |
|---|---|---|---|---|---|
| buried, non-titratable | 23 | +0.18 | 0.19 | 1.06 | 0.00 |
| buried, titratable | 21 | +0.10 | 0.85 | 5.13 | 0.14 |
| exposed, non-titratable | 22 | +0.06 | 0.11 | 0.26 | 0.00 |
| exposed, titratable | 21 | −0.11 | 0.13 | 0.50 | 0.00 |

**Across all 43 exposed-site mutations the mean pH 2.5 → 5.0 shift is 0.12 kcal/mol, the largest
single shift is 0.50 kcal/mol, and none of the 43 exceeds FoldX's own RMSE.** Restricting to
exposed *titratable* residues — the Asp/Glu/Lys/Arg/Tyr the 2025 revision specifically added — makes
no difference: 0.13 kcal/mol.

For scale: the effect the project proposed to build a figure and a hypothesis on is roughly **one
tenth** of the tool's point RMSE (1.25 kcal/mol) and **one twenty-fourth** of its 95 % prediction
interval (2.9 kcal/mol).

The pH response is not absent everywhere — it concentrates exactly where the model has something to
work with, buried titratable residues, where burial raises the ionisation cost. Asp168 (rel. SASA
0.00) moves 0.94–1.71 kcal/mol across D→A/E/K/R/N, and the catalytic Glu188 moves 0.19–0.75. But
those are buried residues, and Fukushima 2017 places the convergent substitutions at the surface.
Where the project needs signal, there is none to separate from noise.

This is an independent, empirical confirmation of the documentary finding in §3.2: the pH module
did not reach statistical significance on its developers' own validation set, and it was validated
on a single histidine. Running it at three pH values and differencing produces numbers, but not
evidence.

---

## 2. The premise problem: Fukushima 2017 already answered this

The project's justification chain is: Fukushima et al. 2017 found convergent substitutions →
their structural/functional consequences are unknown → CarnivorEnzyme measures those consequences
with ΔΔG. But the 2017 paper does report where those substitutions sit and what it takes that to
mean. Verbatim, from the main text:

> "Convergent positions **do not overlap with or cluster around catalytically essential amino
> acids** (Supplementary Fig. 8). Instead, they tend to be located at **exposed positions** to an
> extent comparable to divergent substitutions (Fig. 3c), despite the prediction that more exposed
> positions result in lower convergence probability. Exposed sites are structurally less
> constrained, and substitutions in such sites are likely to **change their interactions with other
> molecules in solution, rather than changing protein conformation**."

Three consequences follow, and each cuts a planned analysis:

1. **ΔΔG of folding measures conformational/thermodynamic consequence.** The source paper states
   these substitutions act on solution interactions instead. §1 confirms the tool has almost no
   dynamic range there. The stability axis is measuring the thing the paper says is not happening.
2. **Catalytic-residue pKa is the wrong target too.** Convergent positions do not cluster near
   catalytic residues, so a CpHMD campaign aimed at the catalytic Asp/Glu dyad is not measuring the
   convergent substitutions' effect — it is measuring a different set of residues.
3. **"Interactions with other molecules in solution" is a testable claim, and it maps onto surface
   electrostatics** — Phase 6, currently a stub. See §5.

Note also that Fukushima 2017 found significant convergence in **GH19 chitinases, purple acid
phosphatases and RNase T2** — not aspartic proteases — and that the convergent branch pairs are
predominantly *Cephalotus* ↔ *N. alata*. That is two origins for most of the signal.

---

## 3. Citation integrity — this doc set contains fabricated references

Both verification passes independently found the same signature: **the DOI or arXiv ID is correct
while the author name is invented.** That is what citations reconstructed from search snippets look
like, not citations that were read. Five confirmed cases:

| Cited as | Reality |
|---|---|
| "Botte et al. (2025)" — FoldX revision (CLAUDE.md, README, both archived FOLDX docs) | **No author named Botte.** Delgado J, Reche R, Cianferoni D, Orlando G, van der Kant R, Rousseau F, Schymkowitz J, Serrano L (2025) "FoldX force field revisited, an improved version," *Bioinformatics* 41(2):btaf064. DOI correct. |
| "Cai et al. (2025) SPURS: Stability prediction using rewiring" | **No author named Cai; title fabricated.** |
| "Cao et al. (2025) … *Nat. Commun.* **16**" | **No author named Cao; wrong volume; wrong title.** DOI, bioRxiv ID and GitHub URL in that entry are correct. |
| "Ressa et al. (2026)" — JanusDDG | **No author named Ressa.** Barducci G, Rossi I, Codicé F, Rollo C, Repetto V, Pancotti C, Iannibelli V, Sanavia T, Fariselli P (2026) *Commun Biol* 9(1):494. |
| "Gapsys et al. 2022, *JCTC* 18:6320" **and** "Guo et al. 2022, *JCTC* 18:6320" — GROMACS CpHMD | **Both fabricated, sharing one bogus locator.** JCTC 18:6320 is Dey, MacAinsh & Zhou on intrinsically disordered protein backbone dynamics. Real refs in §8. |

Two further defects, not fabrications but wrong:

- **PypKa DOI `10.1093/nar/gkae412`** resolves to Zhou et al., *DDMut-PPI* — a different paper.
  Correct is `gkae255`.
- **"CatOpt — Wang et al. 2025"**: first author is Qiu S. Journal, DOI, R² = 0.479 and the ~7 %
  activity gain are all correct.

The one citation the August-2025 pass got **right** and the May docs got wrong is SPURS:
**Li Z, Luo Y (2025), *Nat Commun* 17:891** is correct. (Note it has since carried an Author
Correction, *Nat Commun* 17:2448, 2026-03-16.)

**Action:** every remaining citation in `CLAUDE.md`, `README.md` and `literature_review/` that has
not been resolved against Crossref should be treated as unverified until it is. This is TODO Tier 1.

### 3.1 The FoldX benchmark numbers in CLAUDE.md are real but mislabelled

CLAUDE.md's Method Justification Table says the 2025 revision improved "Pearson R from 0.693 →
0.711, RMSD to 1.238 kcal/mol". Those numbers appear in the paper, but they describe **v10M — the
median of five FoldX runs** (`numberOfRuns=5`), not the force-field revision. The paper's own
headline for the revision is **R 0.693 → 0.706, RMSE 1.252 kcal/mol** on the 2,275-mutation VFSD.
AUC 0.78 is correct.

Separately, `archive/FOLDX_REVIEW_2026-05-12.md`'s claim that "MAVISp's March 2026 evaluation gives
FoldX PCC ~0.40 on Megascale" is wrong three ways: no March-2026 MAVISp release exists, MAVISp has
never benchmarked FoldX against Tsuboyama (it classifies rather than correlating), and 0.40 appears
in no source. The real FoldX-vs-Megascale numbers are from Elwes et al. (bioRxiv, v2 2026-03-05):
**r ≈ 0.30–0.31 per protein, ρ = 0.45 pooled across 58 proteins, rising to r ≈ 0.75** only after
median-aggregation across structures and outlier removal.

### 3.2 The FoldX pH module is far less validated than the docs imply

From the Delgado et al. 2025 paper itself:

- pH dependency for Asp, Glu, Lys, Arg, Tyr is **real** (previously only His and Cys), and is
  implemented as intrinsic experimental pKa + a pairwise-electrostatic environmental correction +
  a burial-based ionisation cost. So it *is* structure-specific — the August-2025 assessment was
  right about that. It is **not** Poisson–Boltzmann, has no Monte Carlo over protonation
  microstates, and the conformation does not respond to protonation.
- It was introduced in **v6, bundled with pi–pi stacking**, and the paper states:
  **"The Wilcoxon signed-rank test of v5 and v6 for the VFSD is not statistically significant
  (P = .373)."** The pH increment did not reach significance on the developers' own validation set.
- The **entire** pH validation is Figure 3: the pH dependence of ΔΔG for **barnase His18→Ala** — one
  mutation, one protein, a *histidine* (pKa ≈ 6.5). There is no acidic-pH stability validation set
  and **no published valid pH range**.

Reporting three ΔΔG values per site at pH 2.5 / 3.5 / 5.0 manufactures three numbers where the tool
supports one, and the differences between them are unvalidated in exactly the regime the project
cares about.

### 3.3 CpHMD cannot carry the pH claim as currently specified

CLAUDE.md describes Phase 5E as "GROMACS 2024.2 + phbuilder | Native λ-dynamics CpHMD …
conda-installable." All three parts are wrong:

- **Constant-pH MD is in no released GROMACS**, through 2026.1. It exists as a fork of **GROMACS
  2021** at `gitlab.com/gromacs-constantph/constantph`, whose README calls the code "preliminary…
  not of production quality" and asks researchers to **avoid publishing results obtained with the
  current version**. Jansen et al. 2024 still call it a beta. Force field is locked to a modified
  CHARMM36m; there is no GPU update support in `gmx mdrun`.
- **Three pH points cannot yield a pKa.** The GROMACS CpHMD authors' own titrations used 5 replicas
  × 50 ns × **15 pH values** (cardiotoxin V) and 5 × 75 ns × **21 pH values** (HEWL). `config.yaml`
  budgets 3 replicas × 50 ns × **3 pH values**. That is protonation fractions at three points, not
  a titration curve — which is in fact what `parse_cphmd.py` outputs (`protonated_fraction`), so
  the script is honest and the plan around it is not.
- **Cost is underscoped 8–17×.** A real titration on a ~400-residue enzyme is ~580–1,220 GPU-hours
  per protein, against the ~64–100 the current budget implies.
- **Accuracy on the directly relevant case is poor.** Hofer et al. 2020 report CpHMD on aspartic
  protease dyads: pepsin good (2.1/5.5 vs. experimental 1.57/5.02), but **HIV-1 protease's lower
  dyad pKa off by 3.6–4.2 pH units**. That error exceeds the entire pH 1.9–6.0 range of interest.

This does not mean CpHMD is worthless — it means it cannot be *the* pH axis, and the 2026-08-25
argument that dropping FoldX was safe "because CpHMD covers pH" does not hold. It also rests on a
category error: `parse_cphmd.py` emits protonated fractions and pKa; FoldX emits folding free
energy. Neither substitutes for the other.

---

## 4. Predictor-by-predictor assessment

### 4.1 The three that survive

**FoldX 5.1 — primary.** The decisive evidence is Pancotti et al. 2022 (*Brief Bioinform*
23(2):bbab555), the S669 benchmark across 21 predictors with an RSA split: *"Most predictors (even
sequence-based) show much lower Pearson correlations on surface residues, **with the exception of
FoldX**, and to a lower extent PremPS and INPS3D."* Given that this project's sites are surface
sites, that exception is the whole argument. Its known weakness is the mirror image: Rollo et al.
2023 (*Genes* 14(12):2228) found FoldX shows "a discernible decrease" in correlation on *modelled*
backbones, where ACDC-NN and DDGun are most robust — and this pipeline feeds it AF3 predictions.
That tension is real and is the reason for a cross-check rather than a sole axis.

**RaSP — the one cross-check.** Blaabjerg et al. 2023 (*eLife* 12:e82593) report crystal-to-crystal
RaSP concordance ρ̄ = 0.97 versus **crystal-to-AlphaFold ρ̄ = 0.94**. That is the strongest published
AlphaFold-robustness evidence for any ΔΔG tool, and it is precisely the axis on which FoldX is
weakest. Cheap, fast, gives saturation maps. *Caveat: the repo carries no LICENSE file — clear
before redistribution.*

**SPURS — optional third, demoted.** Li & Luo 2025 (*Nat Commun* 17:891) is real, peer-reviewed,
MIT-licensed, and genuinely strong: Megascale Spearman median 0.83 (ThermoMPNN 0.77), Human
Domainome 0.54 (0.49), and the "better on 7 of 8 test sets" claim is verbatim in the paper. But the
reason it was made *sole* axis does not hold: the paper uses AlphaFold structures only **as a
fallback when no PDB ID exists** and reports **no experimental-vs-AF ablation**. "AF-native by
design" overstates it. It has no pH input, is trained on ≤72-residue domains, and has 15 citations.

### 4.2 The ones to cut, and why

- **CatOpt and EpHod.** Both are *whole-sequence* pH-optimum regressors — sequence in, one scalar
  out. Mutation effects are obtained by re-running on the mutant and differencing. CatOpt's RMSE is
  0.833 pH units, EpHod's 1.25. Differencing two such predictions carries that error on each side,
  so a single convergent substitution cannot produce a signal above noise. EpHod's own paper warns
  the approach "may not extend to predicting effects of mutations on proteins." Using either as a
  per-site scorer is a methodological error, independent of the citation errors.
- **PypKa in place of PROPKA.** On the pKAI benchmark (736 experimental values, PKAD): PROPKA RMSE
  1.11, PypKa 1.07, and **assuming the solution pKa and doing nothing at all: 1.09**. The swap buys
  0.04 pH units. Not worth the integration work. Keep PROPKA via PDB2PQR, and state its real error —
  0.79 for Asp/Glu overall (Olsson et al. 2011), but **1.37 for buried carboxylates** (Chen et al.
  2022), which is where a catalytic dyad sits.
- **ThermoMPNN, PROSTATA, Pythia, DDMut, JanusDDG, DDGun3D, Stability Oracle** as additional axes:
  see §4.3.
- **Nothing from 2026 as "state of the art."** The 2026 crop is thin and uncited; the two
  best-cited new entries (CATH-ddG, USP-ddG) address protein–protein binding, not folding.
  GraphESMStable and ProStab-Former are real but are unreviewed, zero-citation, single-institution
  preprints — do not cite them as methods.

### 4.3 Why the "4-axis stability ensemble" is rejected

The archived May-2026 docs proposed FoldX × SPURS × RaSP × CatOpt with DDGun3D as a baseline, on
the reasoning that a site flagged by all axes is "very high-confidence." That reasoning does not
survive contact with how these predictors were built.

They are trained on overlapping data (ProTherm, S2648, VariBench, Tsuboyama Megascale) and share a
documented directional bias. Usmanova et al. 2018 (*Bioinformatics* 34:3653) measured it by the
self-consistency test — predicting A→B and B→A on both structures, which must sum to zero:
FoldX **+0.74 ± 0.05 kcal/mol** bias, Eris +1.25, Rosetta +2.08, I-Mutant +0.80. The cause is
explicit: *"programs are trained on experimental datasets which have much more destabilizing
mutations than stabilizing ones… the algorithm will predict more destabilizing mutations,
reflecting the tendencies of the training dataset."* Pucci et al. 2018 found the same in the
direct/inverse framing; Fang 2020 found >70 % sign-inconsistency across five ML predictors;
Savojardo et al. 2021 showed the fix is antisymmetry by construction, not consensus.

**Agreement between predictors that share training data and share a bias measures the shared bias,
not the truth.** Stacking four of them produces a more confident wrong answer, not a more reliable
one. The instinct to prefer one well-chosen axis over an ensemble was correct; what was wrong was
which axis got chosen and how its output would be used.

One further calibration: Montanucci et al. 2019 (*Bioinformatics* 35:1513) derive a **natural upper
bound of Pearson 0.7–0.8** on this task given experimental uncertainty in the reference data, and
note that *"higher Pearson correlations might be indicative of overtraining."* Any predictor
advertising better than that on a standard set deserves suspicion rather than adoption.

---

## 5. Where the pH claim should actually live

The pH question is the scientific heart of the project and it does not belong on the stability
axis. Three independent lines converge on surface electrostatics:

1. **Fukushima 2017's own interpretation** — convergent substitutions change "interactions with
   other molecules in solution" at exposed sites.
2. **The acidophile-adaptation literature.** Fushinobu et al. 1998 (*Protein Eng* 11:1121) on
   *Aspergillus kawachii* xylanase C (pH optimum 2.0): the adaptation is "numerous acidic residues
   concentrated on the surface" plus one 2.8 Å hydrogen bond to the catalytic Glu — and **D37N
   alone shifts the pH optimum from 2.0 to 5.0**. Joshi et al. 2000 (*JMB* 299:255) show the same
   mechanism in *B. circulans* xylanase via reverse protonation. Ohara et al. 2014 (*JBC*
   289:24499) converted an acidophilic carboxylesterase to alkaliphilic with four mutations,
   explained by "a large area of negative charge… around the active site." None of this is visible
   to ΔΔG of folding.
3. **Talley & Alexov 2010** (*Proteins* 78:2699), the one paper that partly defends the stability
   framing, states the caveat plainly: *"only activity is biologically important, while stability
   may not be crucial for the corresponding reaction,"* and notes cases where the two pH optima
   "differ by several pH units."

**Recommendation:** promote Phase 6 (`run_electrostatics.py`, PDB2PQR + PROPKA + APBS at pH 2.5 vs
5.0) from deprioritised stub to the primary pH axis, framed as a surface-charge / active-site
electrostatic-potential comparison across orthologs. That is a claim the tools can actually support
and that the literature has an established mechanism for.

### 5.1 The phenotype side is underpowered — state it or fix it

`config/species.yaml` carries `pitcher_pH` for **9 of 30 species**, spanning only 2.5–5.5, and
those 9 fall in **5 `carnivory_origin` groups** with four of them inside *Nepenthes*. A PGLS of any
enzyme property against pitcher pH is therefore n ≈ 9 with strong phylogenetic clustering. Either
expand the phenotype table from published sources (Adlassnig et al. 2011 is the obvious candidate)
or report the comparative analysis as descriptive rather than inferential. Do not report a p-value
from this as a headline.

---

## 6. Decisions

1. **Reverse `PROJECT_PLAN.md` §2.2.** FoldX is reinstated as the primary stability predictor;
   SPURS is not the sole axis. The August-2025 reasoning failed on two points: FoldX's pH ability
   is not redundant with CpHMD (different observables), and CpHMD cannot deliver what was assumed.
2. **Re-scope the stability axis from classification to ranking.** Drop the
   `ddg_destabilizing: 1.0` / `ddg_stabilizing: -0.5` kcal/mol thresholds for surface sites; report
   within-family percentile rank and the distribution, never a per-site kcal/mol value. §1 shows the
   thresholds are degenerate at the sites of interest.
3. **Run FoldX at a single reference pH.** Retain the pH sweep only as a clearly-labelled
   supplementary sensitivity check, with §3.2's caveats stated in the Methods.
4. **Add RaSP as the sole cross-check.** Where FoldX and RaSP disagree, report it as uncertainty;
   there is no basis for picking a winner.
5. **Move the pH claim to Phase 6 electrostatics.** Promote it in the build order.
6. **Downgrade CpHMD (Phase 5E) to exploratory**, contingent on the fork's publication caveat being
   resolved and on a real titration budget. It cannot be cited as the project's pH result as built.
7. **Cut CatOpt, EpHod, the PypKa swap, and the multi-predictor ensemble** entirely.
8. **Keep `plddt_exclude: 50`.** The proposed tightening to 70 in the archived May docs was
   justified by a "field standard" that does not exist — the 70/90 bands are AlphaFold DB *display*
   colours (Varadi et al. 2022), never validated against ΔΔG accuracy. The only empirically-grounded
   threshold in a ΔΔG context is Akdel et al. 2022's <50, which is what the repo already has.
9. **Sweep every remaining citation against Crossref** before any of it reaches a manuscript.

---

## 7. What this changes in the repo

| File | Change |
|------|--------|
| `PROJECT_PLAN.md` §2.2 | Rewrite: FoldX primary + RaSP cross-check + SPURS optional, replacing "SPURS only, FoldX dropped entirely" |
| `PROJECT_PLAN.md` §4 | Phase 4 gate rewritten for ranking, not threshold classification; Phase 6 promoted; Phase 5E downgraded |
| `PROJECT_STATUS.md` §2 | Method-stack table updated; `run_foldx_*.py` no longer "slated for removal" |
| `CLAUDE.md` §1 | FoldX row: correct the citation to Delgado et al. 2025, correct the 0.711/1.238 attribution, state the pH-validation limits. CpHMD row: correct to Aho et al. 2022 / Jansen et al. 2024 and remove "native"/"conda-installable" |
| `CLAUDE.md` §7 Phase 4/5E | Rewrite gates per §6 |
| `config/config.yaml` | `foldx.ph_values` → single reference pH (sweep moved to a clearly-named supplementary key); thresholds annotated as buried-site-only; `cphmd.run_on_quadrants` re-pointed off the FoldX×EVmutation labels |
| `TODO.md` Tier 1 | T1.2 rewritten: complete `run_foldx_*.py` rather than delete them; add `score_rasp.py`; add the Crossref citation sweep; add the dN/dS + PGLS gaps carried over from the archived paper-structure doc |
| `archive/` | Six superseded root docs moved here on 2026-09-01 with headers |

---

## 8. Sources

Every entry below was resolved against Crossref, PubMed/PMC, Europe PMC, OpenAlex, bioRxiv/openRxiv
or the GitHub API on 2026-09-01.

**Stability predictors**
- Delgado J, Reche R, Cianferoni D, Orlando G, van der Kant R, Rousseau F, Schymkowitz J, Serrano L (2025) FoldX force field revisited, an improved version. *Bioinformatics* 41(2):btaf064. doi:10.1093/bioinformatics/btaf064
- Li Z, Luo Y (2025) Generalizable and scalable protein stability prediction with rewired protein generative models. *Nat Commun* 17(1):891. doi:10.1038/s41467-025-67609-4 *(Author Correction: Nat Commun 17:2448, 2026)*
- Blaabjerg LM, Kassem MM, Good LL, et al. (2023) Rapid protein stability prediction using deep learning representations. *eLife* 12:e82593. doi:10.7554/eLife.82593
- Dieckhaus H, Brocidiacono M, Randolph NZ, Kuhlman B (2024) Transfer learning to leverage larger datasets for improved prediction of protein stability changes. *PNAS* 121(6). doi:10.1073/pnas.2314853121
- Diaz DJ, Gong C, Ouyang-Zhang J, et al. (2024) Stability Oracle. *Nat Commun* 15. doi:10.1038/s41467-024-49780-2
- Barducci G, Rossi I, Codicé F, et al. (2026) JanusDDG. *Commun Biol* 9(1):494. doi:10.1038/s42003-026-09632-9
- Cagiada M, Ovchinnikov S, Lindorff-Larsen K (2025) *Protein Science* 34(1):e5233. doi:10.1002/pro.5233 *(absolute ΔG, not ΔΔG)*
- Tsuboyama K, Dauparas J, Chen J, et al. (2023) Mega-scale experimental analysis of protein folding stability. *Nature* 620:434–444. doi:10.1038/s41586-023-06328-6

**Benchmarks, bias and limits**
- Pancotti C, Benevenuta S, Birolo G, et al. (2022) *Brief Bioinform* 23(2):bbab555. doi:10.1093/bib/bbab555 — S669; FoldX is the exception at surface residues
- Zheng Z, et al. (2024) *Protein Science* 33(1). doi:10.1002/pro.4861 — S4038; all methods better at buried than exposed
- Usmanova DR, Bogatyreva NS, Ariño Bernad J, et al. (2018) *Bioinformatics* 34(21):3653–3658. doi:10.1093/bioinformatics/bty340 — self-consistency bias; FoldX +0.74 kcal/mol
- Pucci F, Bernaerts KV, Kwasigroch JM, Rooman M (2018) *Bioinformatics* 34(21):3659–3665. doi:10.1093/bioinformatics/bty348
- Fang J (2020) *Brief Bioinform* 21(4):1285–1292. doi:10.1093/bib/bbz071; and the rebuttal, Savojardo C, Martelli PL, Casadio R, Fariselli P (2021) *Brief Bioinform* 22(1):601–603. doi:10.1093/bib/bbz168
- Montanucci L, Martelli PL, Ben-Tal N, Fariselli P (2019) A natural upper bound to the accuracy of predicting protein stability changes upon mutations. *Bioinformatics* 35(9):1513–1517. doi:10.1093/bioinformatics/bty880
- Sapozhnikov Y, Patel JS, Ytreberg FM, Miller CR (2023) *BMC Bioinformatics* 24:426. doi:10.1186/s12859-023-05537-0 — ±2.9 kcal/mol 95% prediction interval
- Rollo C, et al. (2023) *Genes* 14(12):2228. doi:10.3390/genes14122228 — predictor sensitivity to modelled structures
- Broom A, Jacobi Z, Trainor K, Meiering EM (2017) *JBC* 292(35):14349–14361. doi:10.1074/jbc.M117.784165
- Elwes C, Alcraft R, Lister H, Smith PA, Shorthouse D, Hall BA (2026) Validating folding energy estimates as a method for variant interpretation. bioRxiv doi:10.1101/2025.11.09.687451 *(preprint)*

**AlphaFold as ΔΔG input**
- Akdel M, Pires DEV, Porta Pardo E, et al. (2022) A structural biology community assessment of AlphaFold2 applications. *Nat Struct Mol Biol* 29(11):1056–1067. doi:10.1038/s41594-022-00849-w
- Pak MA, Markhieva KA, Novikova MS, et al. (2023) *PLOS ONE* 18(3):e0282689. doi:10.1371/journal.pone.0282689 — ΔpLDDT vs ΔΔG r = −0.17
- Varadi M, Anyango S, Deshpande M, et al. (2022) *NAR* 50(D1):D439–D444. doi:10.1093/nar/gkab1061 — the pLDDT bands are display categories

**pH, pKa and constant-pH MD**
- Aho N, Buslaev P, Jansen A, Bauer P, Groenhof G, Hess B (2022) Scalable constant pH molecular dynamics in GROMACS. *JCTC* 18(10):6148–6160. doi:10.1021/acs.jctc.2c00516
- Buslaev P, Aho N, Jansen A, Bauer P, Hess B, Groenhof G (2022) Best practices in constant pH MD simulations. *JCTC* 18(10):6134–6147. doi:10.1021/acs.jctc.2c00517
- Jansen A, Aho N, Groenhof G, Buslaev P, Hess B (2024) phbuilder. *JCIM* 64(3):567–574. doi:10.1021/acs.jcim.3c01313
- Hofer F, Kraml J, Kahler U, Kamenik AS, Liedl KR (2020) *JCIM* 60(6):3030–3042. doi:10.1021/acs.jcim.0c00190 — CpHMD on aspartic protease dyads
- Olsson MHM, Søndergaard CR, Rostkowski M, Jensen JH (2011) PROPKA3. *JCTC* 7(2):525–537. doi:10.1021/ct100578z
- Chen AY, Lee J, Damjanovic A, Brooks BR (2022) *JCTC* 18(4):2673–2686. doi:10.1021/acs.jctc.1c01257 — buried-carboxylate pKa error
- Reis PBPS, Clevert D-A, Machuqueiro M (2024) PypKa server. *NAR* 52(W1):W294–W298. doi:10.1093/nar/gkae255
- Qiu S, Wang N-K, Lu Y, Gong J-S, Shi J-S, Yang A (2025) CatOpt. *ACS Synth Biol* 14(12):4897–4906. doi:10.1021/acssynbio.5c00679
- Gado JE, Knotts M, Shaw AY, et al. (2025) Machine learning prediction of enzyme optimum pH. *Nat Mach Intell* 7(5):716–729. doi:10.1038/s42256-025-01026-6

**Acid adaptation mechanism**
- Fushinobu S, Ito K, Konno M, Wakagi T, Matsuzawa H (1998) *Protein Eng* 11(12):1121–1128. doi:10.1093/protein/11.12.1121
- Joshi MD, Sidhu G, Pot I, Brayer GD, Withers SG, McIntosh LP (2000) *JMB* 299:255–279. doi:10.1006/jmbi.2000.3722
- Ohara K, Unno H, Oshima Y, et al. (2014) *JBC* 289(35):24499–24510. doi:10.1074/jbc.M113.521856
- Lin Y, Fusek M, Lin X, Hartsuck JA, Kezdy FJ, Tang J (1992) *JBC* 267:18413–18418. PMID 1526982
- Talley K, Alexov E (2010) *Proteins* 78(12):2699–2706. doi:10.1002/prot.22786

**Project foundation**
- Fukushima K, Fang X, Alvarez-Ponce D, et al. (2017) *Nat Ecol Evol* 1:0059. doi:10.1038/s41559-016-0059 — §2 quotation is from the main text
- Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013) *PLOS ONE* 8:e80635 — max ASA normalisation used in §1

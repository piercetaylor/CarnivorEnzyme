# CarnivorEnzyme — Project Plan

> Date: 2026-08-26 (§2.2 and §4 revised 2026-09-01 — see `audit/15_stability_predictor_audit.md`)
> Supersedes: `PROJECT_RESTRUCTURE_2026-05-12.md` (archived, see `archive/PROJECT_RESTRUCTURE_2026-05-12.md`)
> Companion documents: [PROJECT_STATUS.md](PROJECT_STATUS.md), [TODO.md](TODO.md)
>
> This is a **decision record**, not a task list — it captures the scientific framing and the
> method-stack choices this project has committed to, with the reasoning behind each. For what's
> actually built vs. stubbed, see PROJECT_STATUS.md. For the actionable, gated backlog, see
> TODO.md.

---

## 1. The scientific question

Carnivory evolved independently in multiple angiosperm lineages (see §3 below for the corrected
count). Each lineage co-opted pre-existing hydrolase gene families — chitinases, acid
phosphatases, RNases, aspartic proteases — for extracellular digestion. Fukushima et al. (2017,
*Nat. Ecol. Evol.* 1:0059) established that the same amino acid substitutions recur at homologous
alignment positions across phylogenetically independent lineages at a rate exceeding neutral
expectation. **The structural and functional consequences of these convergent substitutions remain
uncharacterized** — this is the question CarnivorEnzyme answers, and it's explicitly flagged as
open by the 2025 Albert/Fukushima review (see CLAUDE.md §0).

Full scientific framing, target enzyme families, and lineage table live in CLAUDE.md §§0, 3, 4 —
those sections are current (corrected 2026-08-22 for MEROPS codes, GH19 class split, and the
carnivory-origin count) and are not re-derived here. This document adds the method-stack decisions
CLAUDE.md's tool table has not yet caught up to.

---

## 2. Method-stack decisions (resolved 2026-08-25; §2.2 reversed 2026-09-01)

CLAUDE.md's own addendum instructs: "This project must not become a confirmation loop... every
method choice must be challenged against alternatives." The addendum already did this once (FoldX
alone → FoldX + EVmutation). A second round was needed because `PROJECT_RESTRUCTURE_2026-05-12.md`
had proposed a further pivot (ProSST/SaProt, SPURS/RaSP, EVE) that sat unresolved for over three
months — TASKS.md itself flagged it "STALE, unexecuted." The decisions below close that out, based
on live research (real web search, not memory) conducted 2026-08-25.

### 2.1 Outcome track: computational-only (Outcome A)

No wet-lab collaboration outreach. Target venue tier: MBE/GBE, fully computational. This resolves
the fastest path to a submittable preprint and removes the external-dependency risk a wet-lab
partnership would add. (Old TASKS.md T4.1 — wet-lab outreach — is dropped, not carried forward.)

### 2.2 Stability ΔΔG: FoldX primary + RaSP cross-check (REVISED 2026-09-01)

> **This section was rewritten on 2026-09-01.** The prior decision — "SPURS only, FoldX dropped
> entirely" — is **reversed**. See `audit/15_stability_predictor_audit.md` for the full evidence,
> including first-hand FoldX 5.1 runs on the 7ZVA crystal and two independent citation-verification
> passes. The superseded reasoning is preserved at the end of this section.

**Decision:** **FoldX 5.1** is the primary stability predictor, **RaSP** is the single cross-check,
and **SPURS** is optional and never sole. Stability is reported as a within-family **rank**, not as
per-site kcal/mol against fixed thresholds. The pH sweep is demoted from a headline result to a
labelled supplementary sensitivity check, and the **pH claim moves to Phase 6 electrostatics**.

**Why FoldX comes back:**

- **It is the one predictor that holds up at surface residues**, which is where this project's
  substitutions are. Pancotti et al. 2022 (*Brief Bioinform* 23(2):bbab555, the S669 benchmark
  across 21 predictors with an RSA split): "Most predictors (even sequence-based) show much lower
  Pearson correlations on surface residues, **with the exception of FoldX**." Zheng et al. 2024
  (*Protein Sci* 33(1)) confirm every method is worse at exposed than buried sites.
- **Its counter-argument is real and is why RaSP joins it.** Rollo et al. 2023 (*Genes* 14:2228)
  found FoldX degrades measurably on *modelled* backbones — and this pipeline feeds it AF3
  predictions. RaSP is the only tool with published crystal-vs-AlphaFold concordance
  (ρ̄ = 0.94, against 0.97 crystal-to-crystal; Blaabjerg et al. 2023). Where the two disagree,
  report uncertainty; there is no basis to pick a winner.
- **It is already licensed and installed** (FoldX 5.1, verified 2026-09-01), so the reinstatement
  costs integration work only.

**Why SPURS is demoted, not dropped:** SPURS is real, peer-reviewed and strong (Li & Luo 2025,
*Nat Commun* 17:891 — that citation was correct; the May docs' "Cao"/"Cai" attributions are
fabricated author names). But the reason it was made *sole* axis does not survive checking: the
paper uses AlphaFold structures **only as a fallback where no PDB ID exists** and reports **no
experimental-vs-AF ablation**, so "AF-native by design" overstates it. It also has no pH input, is
trained on ≤72-residue domains, and has 15 citations.

**Why the ensemble stays rejected.** The instinct that one well-chosen axis beats a stack was
right, and it still holds — but for a sharper reason than recency. These predictors share training
data (ProTherm, S2648, VariBench, Megascale) and a documented directional bias toward destabilizing
predictions (Usmanova et al. 2018 measure FoldX at **+0.74 kcal/mol** by the self-consistency test;
Pucci et al. 2018 and Fang 2020 find the same). Agreement between them measures the shared bias,
not the truth. RaSP is added as a cross-check on *structure sensitivity*, a failure mode FoldX has
and RaSP demonstrably does not — not as a second vote.

**What the ΔΔG axis may and may not claim.** Direct measurement on 7ZVA (87 mutations, buried and
exposed, FoldX 5.1) gives mean |ΔΔG| of **6.1–7.1 kcal/mol at buried positions but 0.39–0.42 at
exposed ones**, with 0 of 43 exposed mutations exceeding the published ±2.9 kcal/mol 95 % prediction
interval and 2 of 43 exceeding the 1.25 kcal/mol RMSE. Against `config.yaml`'s
`ddg_destabilizing: 1.0`, every buried site classifies as destabilizing and essentially every
surface site as neutral — the axis is degenerate at the sites of interest. Hence: **ranked output,
distributions across families, no per-site kcal/mol claims at surface positions.**

**Why the pH sweep is demoted.** Two independent reasons, one documentary and one measured:

- The FoldX pH module was introduced in v6 **bundled with pi–pi stacking**, and the paper states
  "The Wilcoxon signed-rank test of v5 and v6 for the VFSD is not statistically significant
  (P = .373)." Its entire validation is one figure: barnase His18→Ala — one mutation, one protein,
  a *histidine*. There is no acidic-pH validation set and no published valid pH range.
- Measured on 7ZVA, the **pH 2.5 → 5.0 shift at exposed sites averages 0.12 kcal/mol** (max 0.50);
  none of 43 mutations clears the tool's own RMSE. The planned Fig. 3 / H1 effect is ~10× smaller
  than the error bar.

**Why CpHMD cannot absorb the pH claim** (the load-bearing assumption of the superseded decision):
constant-pH MD is in **no released GROMACS** through 2026.1 — it is a fork of GROMACS 2021 whose
README asks researchers to avoid publishing results from it; three pH points cannot produce a pKa
(the developers' own titrations used 15–21); the budget is underscoped 8–17×; and on the directly
relevant case, Hofer et al. 2020 missed HIV-1 protease's lower catalytic-dyad pKa by 3.6–4.2 pH
units. Beyond feasibility, the substitution was a category error: `parse_cphmd.py` emits protonated
fractions and pKa, FoldX emits folding free energy. Neither replaces the other.

**Where the pH claim goes instead:** Phase 6 electrostatics (PDB2PQR + PROPKA + APBS at pH 2.5 vs
5.0), promoted from deprioritised stub to primary pH axis. This is what Fukushima 2017's own
interpretation points to — convergent positions are surface-exposed, do not cluster near catalytic
residues, and are read as changing "interactions with other molecules in solution, rather than
changing protein conformation" — and it matches the published mechanism of acid adaptation
(Fushinobu et al. 1998; Ohara et al. 2014: acidic residues concentrated on the surface, negative
potential around the active site).

**Consequence:** the old Phase 4 (`run_foldx_repair.py`, `run_foldx_scan.py`, `parse_foldx.py`,
`foldx.smk`) is **completed, not deleted**. `config.yaml` keeps a `foldx:` block with a single
reference pH; `score_rasp.py` is added; a `spurs:` block is optional. CLAUDE.md's Method
Justification Table keeps FoldX with a corrected citation and honest limitations.

<details>
<summary>Superseded reasoning (2026-08-25) — preserved for the record</summary>

The 2026-08-25 decision made SPURS the sole stability predictor and removed FoldX entirely, on two
grounds: (a) SPURS is the newest peer-reviewed AF-native method, and (b) FoldX's pH-dependent ΔΔG
is redundant with CpHMD, which is "more rigorous." It also cited FoldX at "PCC ~0.40 on Megascale
per MAVISp's March 2026 evaluation."

All three load-bearing claims failed verification. (a) SPURS is not AF-native in the sense claimed —
AF structures are a fallback, with no ablation. (b) CpHMD and FoldX compute different observables,
and the GROMACS CpHMD implementation cannot deliver a publishable pKa as configured. (c) There is no
March-2026 MAVISp release, MAVISp has never benchmarked FoldX against Tsuboyama, and 0.40 appears in
no source; the real figures are r ≈ 0.30–0.31 per protein and ρ = 0.45 pooled (Elwes et al. 2026).

The decision was made in good faith on documents that turned out to contain fabricated citations —
the FoldX revision was attributed throughout to a nonexistent "Botte et al." (real: Delgado et al.
2025), and both GROMACS CpHMD references were invented.
</details>

### 2.3 Evolutionary fitness: ProSST (primary) + SaProt (secondary), ESM-2 demoted to baseline

**Decision:** ProSST becomes the primary sequence/structure-aware evolutionary-fitness scorer,
SaProt secondary (mirroring I/O contract), ESM-2 retained only as a Tier-3 sanity check —
`score_esm2.py` is not deleted, just no longer load-bearing for the quadrant classification.

**Why:** ProSST tops the ProteinGym leaderboard as of the May 2026 review and is structure-aware,
which ESM-2 (sequence-only) is not. **Real cost:** ProSST/SaProt need Foldseek 3Di tokens extracted
from a predicted structure, which creates a hard dependency on Phase 3 (structure prediction)
landing before this axis can run at all — currently Phase 3 is a 100% stub. This is reflected in
the build-order in §4 below: Phase 3 must land before Phase 3A's ProSST/SaProt scoring can run for
real, even though the *rule file* for Phase 3A already exists.

### 2.4 Evolutionary coupling (DCA): EVE, not EVmutation/EVcouplings

**Decision:** EVE (Frazer et al. 2021) replaces EVmutation/EVcouplings as the deep-MSA
evolutionary-constraint axis for families with Neff ≥ 500.

**Why:** better-calibrated per the May 2026 review; EVmutation has more published precedent but EVE's
benchmark numbers are stronger. **Consequence:** `run_evcouplings.py`, `run_evmutation.py`,
`compare_foldx_evmutation.py`, and `evmutation.smk` currently target the wrong tool — TODO.md Tier 1
renames/rebuilds them for EVE rather than completing them as-is. The Neff ≥ 500 gate itself
carries over unchanged (same MSA-depth requirement applies to EVE).

### 2.5 Quadrant classification, re-derived on the new axes

The primary novel analytical output (see §5) is **FoldX ΔΔG × (ProSST/SaProt LLR, informed by EVE
ΔΔE for deep-MSA families)**, with RaSP shown alongside FoldX as a structure-sensitivity
cross-check. `compare_foldx_evmutation.py` is renamed `compare_stability_evolution.py`.

**The stability axis is a within-family percentile rank, not thresholded kcal/mol** (revised
2026-09-01). The old quadrant assumed a meaningful destabilizing/neutral/stabilizing split at
1.0 / −0.5 kcal/mol. Direct measurement shows that split is degenerate at surface positions, where
this project's substitutions are: all 43 exposed-site mutations tested on 7ZVA fall within
±1.6 kcal/mol, mean |ΔΔG| 0.40, against a tool RMSE of 1.25 and a 95 % prediction interval of ±2.9
(`audit/15_stability_predictor_audit.md` §1). Thresholding would label essentially every convergent
site "neutral" and every buried control "destabilizing", which is an artefact of burial, not a
finding about convergence.

Ranking sidesteps this: it asks "is this substitution unusually costly *relative to other
substitutions in this family*", which the tool can support, rather than "does it cross an absolute
kcal/mol boundary", which it cannot at these sites. Quadrant boundaries are therefore set at
within-family percentiles, fixed before looking at which sites are convergent — TODO.md Tier 1
carries the concrete sub-task.

---

## 3. Target families and lineages

Unchanged from CLAUDE.md §§3–4 (already corrected 2026-08-22): 4 Tier-1 families
(`chitinases_gh19_class_iv`, `purple_acid_phosphatase`, `rnase_t2`,
`aspartic_proteases_a1b_homology`), 3 methods-benchmark families (`chitinases_gh19_class_i`,
`nepenthesins`, `neprosins`), single Caryophyllales carnivory origin. Not re-derived here — see
CLAUDE.md directly, it's current.

**Note:** `STATUS_2026-05-12.md` (archived) described a 5-family, 9-origin framing predating the
Aug-22 corrections. That framing is superseded independent of the method-stack pivot — don't revive
it.

---

## 4. Build order and test gates (updated for the new stack)

Phase numbering follows CLAUDE.md §7, with the tool substitutions from §2 above:

| Phase | What | Test gate (updated) |
|-------|------|----------------------|
| 1 | Sequence retrieval | Unchanged — nepenthesin orthogroup must contain the expected accessions and topology (CLAUDE.md Test Gate 1) |
| 2 | Convergence detection | Unchanged — ≥80% recall of Fukushima 2017 Fig. 3a GH19 positions (CLAUDE.md Test Gate 2). **Blocked on TODO.md Tier 0 HPC verification before this gate can be trusted for real data.** |
| 3A | Evolutionary fitness triage | **New gate:** ProSST LLR positive for ≥4/6 known Fukushima 2017 GH19 convergent positions (same bar the old ESM-2-only gate used). Hard-blocked on Phase 3 landing (needs Foldseek 3Di tokens from a predicted structure). |
| 3 | Structure prediction | Unchanged — AF3 vs. 7ZVA TM-score >0.90 (Chai-1 fallback >0.85), CLAUDE.md Test Gate 3 |
| 4A | Ancestral structure | Unchanged — MRCA <30% identity to any single modern sequence, TM-score >0.70 vs. modern ortholog, CLAUDE.md Test Gate 3B. **Also gated on the `extract_ancestor.py` node-naming audit in TODO.md Tier 0.** |
| 4 | Stability ΔΔG | **Revised 2026-09-01.** FoldX 5.1 primary + RaSP cross-check. Gate: (a) FoldX ΔΔG on the 7ZVA crystal vs. the AF3 neprosin prediction must correlate at r > 0.6 across matched positions — the original Test Gate 4, reinstated; (b) FoldX and RaSP rank-correlate at ρ > 0.5 within a family, with disagreement reported rather than resolved. Output is a within-family **rank**, not thresholded kcal/mol — the 1.0/−0.5 cutoffs are degenerate at surface sites (audit/15 §1) |
| DCA | EVE scoring | **New gate, replaces EVmutation sign check:** negative ΔΔE for known pathogenic/deleterious mutations in a well-studied homolog, Neff ≥ 500 gate unchanged |
| 5 | Docking | Unchanged — PQPQLPYP into neprosin 7ZVC active-site cleft, CLAUDE.md Test Gate 5 |
| 5D | FEP | Unchanged — CLAUDE.md Test Gate 5D |
| 5E | CpHMD | **Downgraded 2026-09-01 to exploratory.** It is *not* the pipeline's pH signal. Blocked on: the GROMACS constant-pH fork's own "do not publish from this version" caveat being resolved, and a titration budget with ≥10 pH points (3 cannot yield a pKa). Until then it produces protonated fractions for a case study, not a headline result (audit/15 §3.3) |
| 6 | Electrostatics | **Promoted 2026-09-01 to the primary pH axis.** PDB2PQR + PROPKA + APBS at pH 2.5 vs. 5.0. Gate: surface-charge / active-site potential differences across orthologs are computable and directionally interpretable for ≥1 Tier-1 family. This carries the acid-adaptation claim that FoldX's pH sweep and CpHMD cannot (audit/15 §5) |
| 7 | Expression | Unchanged |
| 8 | Integration | Fig. 5 (quadrant scatter) becomes **FoldX rank × ProSST/SaProt/EVE**, with RaSP rank shown as the cross-check. `generate_figures.py` must plot percentile rank on the stability axis, not kcal/mol |

---

## 5. Novelty and competition framing (tightened 2026-08-25)

CLAUDE.md's addendum §D currently claims: *"No published study has done this on any system, not
just carnivorous plants."* That's an overstatement — live research (fork agent, real web search,
2026-08-25) found:

- **The method itself is mainstream, not novel.** Crossing a stability-ΔΔG predictor against an
  evolutionary-constraint/conservation predictor to classify variants into a 2D grid is standard
  practice in deep-mutational-scanning and clinical variant-effect-predictor interpretation.
  Gerasimavicius, Livesey & Marsh (*Cell Reports* 2021/2022; *Protein Science* 2023) build exactly
  this kind of stability × conservation classification for general (non-convergent) variants,
  including loss-of-function / gain-of-function / dominant-negative categories. Tokuriki–Tawfik-
  style "stability-function tradeoff" framing is standard in protein engineering.
- **Convergent-evolution precedent is thinner but real in spirit.** Toxin-resistance convergence
  studies (Na,K-ATPase / TTX resistance — Brodie & Feldman lineage, *MBE* 2022, *PNAS* 2009)
  repeatedly document convergent substitutions trading resistance against baseline
  stability/function via structural modeling, but none formalize it as a paired-predictor quadrant.
- **The honest, defensible novelty claim:** no published paper crosses a structure-based stability
  predictor against a sequence-based evolutionary-constraint predictor **specifically on convergent
  substitutions across independent lineages**, in a formal quadrant. The method is well-trodden;
  applying it to convergence data is the actual novel contribution.

CLAUDE.md addendum §D should be revised to this narrower framing when Tier 1 lands (cite
Gerasimavicius/Livesey/Marsh and the toxin-resistance-convergence papers as the two nearest
precedent classes) — see TODO.md Tier 1.

Direct precedent comparisons (Butts 2016, Fukushima 2017, the Rubisco FoldX study, Kim 2025 —
CLAUDE.md addendum §D §1–4) are unaffected by the method-stack pivot and remain accurate as-is.

---

## 6. Scope boundaries under Outcome A

- Tier 2 families (cysteine proteases, TLPs, GH17, LTPs, esterases) stay deferred until Tier 1 is
  validated end-to-end — unchanged from CLAUDE.md §3.
- FEP scope: capped at ≤20 priority sites total, gated on the FoldX×ProSST/EVE quadrant
  (`functional_gain` / `stability_function_tradeoff`) rather than the old FoldX×EVmutation gate.
  The numeric cap (20) carries over from the old T0.3 proposal, but the gating *quadrant labels*
  now come from the new axes — re-verify the cap is still sane once real quadrant data exists
  (TODO.md Tier 1/2).
- No wet-lab section in the eventual manuscript (Outcome A, §2.1).

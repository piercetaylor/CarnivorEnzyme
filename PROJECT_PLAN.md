# CarnivorEnzyme — Project Plan

> Date: 2026-08-26
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

## 2. Method-stack decisions (resolved 2026-08-25)

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

### 2.2 Stability ΔΔG: SPURS only, FoldX dropped entirely

**Decision:** SPURS (Li & Luo, *Nat. Commun.* 17:891, 2025/2026) becomes the sole protein-stability
predictor. FoldX is removed from the pipeline — not run in production, not kept as a benchmark
check, not swept across pH.

**Why not stack FoldX + a second DL predictor (SPURS or RaSP), which was the original
`PROJECT_RESTRUCTURE_2026-05-12.md` proposal?** The user's own instinct was that a single
well-chosen axis beats an ensemble here, which prompted two rounds of research:

- **Round 1 (predictor choice):** compared SPURS, ThermoMPNN, RaSP, CatOpt, ddGUn3D, and two
  January-2026 preprints (GraphESMStable, ProStab-Former) on recency, peer-review status,
  AF-native training, and benchmark correlation. SPURS is the newest *peer-reviewed*, AF-native
  method (Pearson/Spearman correlations of 0.54–0.83 depending on benchmark set) with a real
  citation trail (19 citations as of Aug 2026). GraphESMStable/ProStab-Former report higher raw
  numbers but are unreviewed, zero-citation preprints 4–7 weeks old — not reviewer-defensible yet.
  **Correction to this repo's own prior docs:** `FOLDX_ALTERNATIVES_AF3_2026-05-12.md` and
  `FOLDX_REVIEW_2026-05-12.md` mis-cite SPURS as "Cao et al." — the correct citation is Li & Luo.
- **Round 2 (does dropping FoldX lose something real?):** FoldX's pH-dependent ΔΔG (the 4.1/5.1
  revision, Botte et al. 2025) is genuinely structure-specific — it computes an
  environmentally-corrected pKa per titratable residue from local electrostatics and burial, not a
  generic per-amino-acid formula. But it's answering the same question CpHMD (already implemented,
  Phase 5E) answers with a fundamentally more rigorous method — explicit λ-dynamics MD vs. FoldX's
  static empirical correction. Given FoldX's own general-ΔΔG benchmark is weak and inconsistent
  (PCC ~0.40 on Megascale per MAVISp's March 2026 evaluation, up to ~0.71 depending on dataset) and
  it has no dedicated pH-specific validation set, keeping it *specifically for the pH sweep* adds
  less unique value than it appears. **CpHMD becomes the pipeline's sole pH-dependent stability
  signal.**

**Consequence:** the old Phase 4 (`run_foldx_repair.py`, `run_foldx_scan.py`, `parse_foldx.py`,
`foldx.smk`) is not completed — it's replaced. `config.yaml`'s `foldx:` block is removed; a
`spurs:` block is added. CLAUDE.md's Method Justification Table FoldX row is replaced with SPURS.

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

The primary novel analytical output (see §5) is now **SPURS ΔΔG × (ProSST/SaProt LLR, informed by
EVE ΔΔE for deep-MSA families)**, not the old FoldX × EVmutation pairing. `compare_foldx_evmutation.py`
is renamed `compare_spurs_eve.py`; quadrant thresholds need re-deriving against SPURS's score scale
(SPURS does not report kcal/mol the way FoldX does — TODO.md Tier 1 flags this as a concrete
open sub-task, not something to assume carries over numerically).

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
| 4 | Stability ΔΔG | **New gate, replaces old FoldX Test Gate 4:** SPURS ΔΔG on a held-out structure (e.g. neprosin AF3 prediction) should be directionally consistent with literature-known stabilizing/destabilizing mutations at ≥3 test positions — exact quantitative correlation threshold needs picking once SPURS is actually run once (TODO.md Tier 1 sub-task; SPURS's score scale differs from FoldX's kcal/mol, so the old "r > 0.65" threshold doesn't transfer directly) |
| DCA | EVE scoring | **New gate, replaces EVmutation sign check:** negative ΔΔE for known pathogenic/deleterious mutations in a well-studied homolog, Neff ≥ 500 gate unchanged |
| 5 | Docking | Unchanged — PQPQLPYP into neprosin 7ZVC active-site cleft, CLAUDE.md Test Gate 5 |
| 5D | FEP | Unchanged — CLAUDE.md Test Gate 5D |
| 5E | CpHMD | Unchanged — CLAUDE.md Test Gate 5E. **Now also the pipeline's sole pH-resolved stability signal (see §2.2) — its output carries more analytical weight than before.** |
| 6 | Electrostatics | Unchanged |
| 7 | Expression | Unchanged |
| 8 | Integration | Fig. 5 (quadrant scatter) becomes SPURS × ProSST/SaProt/EVE, not FoldX × EVmutation — figure-generation logic in `generate_figures.py` needs to reflect this when built |

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
- FEP scope: capped at ≤20 priority sites total, gated on the SPURS×ProSST/EVE quadrant
  (`functional_gain` / `stability_function_tradeoff`) rather than the old FoldX×EVmutation gate.
  The numeric cap (20) carries over from the old T0.3 proposal, but the gating *quadrant labels*
  now come from the new axes — re-verify the cap is still sane once real quadrant data exists
  (TODO.md Tier 1/2).
- No wet-lab section in the eventual manuscript (Outcome A, §2.1).

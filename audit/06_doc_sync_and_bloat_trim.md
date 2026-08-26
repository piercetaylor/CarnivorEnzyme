# T6: Documentation sync + bloat trim

**Execution note:** this task was attempted twice by a delegated agent and stalled/crashed both
times on transient infrastructure issues (a 600s no-progress watchdog, then a DNS/API
connectivity error) — neither failure was caused by the task logic itself. The first attempt made
zero edits before dying (verified via `git status` before retry). The second attempt completed
`CLAUDE.md`, `README.md`, and `TASKS.md` correctly before dying mid-task. The orchestrating
session verified what had actually landed, then completed the remaining work directly
(`literature_review/LITERATURE_FINDINGS.md`, `literature_review/MASTER_LITERATURE_REGISTRY.md`,
`ORTHOLOGY_AND_FAMILY_SCOPE_2026-05-12.md`, the bloat-trim decision, and the final consistency
sweep) rather than risk a third stall.

## 1. `CLAUDE.md` — updated (by the delegated agent, verified correct)

- Section 3 (Target Enzyme Families): rewritten to match `config/enzyme_families.yaml`'s actual
  current `tier1:` keys (`chitinases_gh19_class_iv`, `chitinases_gh19_class_i`,
  `purple_acid_phosphatase`, `rnase_t2`, `aspartic_proteases_a1b_homology`) plus a new "Methods
  Benchmark" subsection for `nepenthesins`/`neprosins` (moved out of tier1, single-origin
  families, not cross-lineage convergence targets).
- Section 4 (Carnivorous Plant Lineages): added a correction block citing Fleck & Jobson 2023 and
  Fleischmann et al. 2018 verbatim, matching `config/species.yaml`'s header — Caryophyllales
  carnivory is one shared origin, not three.
- Section 10 citation fixed: "Tiew & Goh 2022... 183:23-38" → "Ting TY, Baharin A, Ramzi AB, Ng
  CL, Goh HH (2022)... 183:23-35" with accurate framing (paper claims "a new" family, not formal
  G3 membership; Oda & Wlodawer 2023 added as the separate structural-similarity citation).

## 2. `README.md` and `TASKS.md` — updated (by the delegated agent, verified correct)

- `README.md`: neprosin citation corrected to Ting et al. 2022 with the same accurate framing;
  origin-count language updated ("one shared origin, not four separate ones" — describing the
  four genera Nepenthes/Drosera/Dionaea/Aldrovanda collapsing to one origin, consistent with but
  phrased differently from the "3 origin values → 1" framing used elsewhere; not a factual error).
  One residual cosmetic issue (markdownlint MD012, a duplicated blank line under "## The
  Question") was found and fixed directly by the orchestrator.
- `TASKS.md`: the stale `carnivory_origin` values for `Sarracenia_purpurea` (was 5), and
  `Pinguicula_gigantea` (was 6) corrected to match the current post-T2 numbering (3 and 4
  respectively; Caryophyllales=1, Oxalidales=2, Ericales=3, Lentibulariaceae=4, Byblidaceae=5),
  with an explanatory note added rather than a silent renumber.

## 3. `literature_review/LITERATURE_FINDINGS.md` — updated (by the orchestrator directly)

- The neprosin bibliography entry (§ Structural Classification) corrected: "Tiew YK, Goh HH
  (2022)... 183:23–38" → "Ting TY, Baharin A, Ramzi AB, Ng CL, Goh HH (2022)... 183:23–35", with
  the same accurate G3-vs-U74 framing added to its Key Findings.
- § 4B (Carnivorous Plant Phylogeny Summary) rewritten: previously listed Nepenthes, Drosera, and
  Dionaea as three separate independent origins within Caryophyllales ("third independent origin
  in Caryophyllales"). Replaced with the corrected single-origin framing, citing Fleck & Jobson
  2023 and Fleischmann et al. 2018, and collapsed the origin count from "5+ lineages" to the
  accurate "5 lineages, 4 orders." The one remaining "independent origin from Caryophyllales" hit
  (describing Cephalotus, which genuinely is a separate/independent origin from Caryophyllales)
  was left as-is — it's correct, not a leftover error.

## 4. Two more stale "Tiew & Goh" citations found outside the original checklist

A full-repo grep (step 5 below) turned up two files not covered by the original task list, both
containing live (uncorrected) citations to the wrong author name:

- **`literature_review/MASTER_LITERATURE_REGISTRY.md`**: "Tiew YK, Goh HH (2022)... 183:23–38...
  doi:10.1016/j.plaphy.2022.04.029" — wrong author, wrong page range, AND a different (also
  wrong) DOI than the one already verified elsewhere in this audit trail
  (10.1016/j.plaphy.2022.04.027, confirmed via Europe PMC REST, PMID 35537348). Corrected all
  three simultaneously.
- **`ORTHOLOGY_AND_FAMILY_SCOPE_2026-05-12.md`**: two live "Tiew & Goh 2022" citations (one
  inline, one in the Sources list) — DOI was already correct here, only the author name and page
  range needed fixing. Corrected both.

## 5. Bloat trim — investigated, nothing deleted

The following six root-level files were flagged by an earlier (less thorough) audit pass as
likely-superseded working notes worth trimming: `STATUS_2026-05-12.md`, `PROGRESS.md`,
`PROJECT_RESTRUCTURE_2026-05-12.md`, `temporary_TODO_42526.md`, `FOLDX_REVIEW_2026-05-12.md`,
`FOLDX_ALTERNATIVES_AF3_2026-05-12.md`. Per this task's own instructions ("read first, then
decide, don't delete blindly"), the orchestrator read all six in full before making any decision.

**Finding: none of these are bloat. All six are substantive, still-relevant technical documents,
not stale/duplicated working notes. None were deleted.**

- **`STATUS_2026-05-12.md`** documents a real, apparently-still-open architectural pivot decision
  (ESM-2 → ProSST + SaProt as the primary evolutionary-constraint scorer; three new genomes;
  FEP scope capped at 20 sites) that does **not** appear to have been executed in the current
  codebase — `workflow/scripts/score_esm2.py` still exists under its original name, and no
  `score_saprot.py`/`score_prosst.py` exists. This is a genuine open PI decision, not resolved
  bloat, and is well outside this remediation pass's scope (taxonomy/contamination/loader-bug
  fixes, not tool-stack architecture).
- **`PROJECT_RESTRUCTURE_2026-05-12.md`** is a detailed, coherent proposal to reorganize the
  Snakefile from method-organized phases to question-organized "paper sections" (A–G plus
  Supplementary), with an explicit "§8 Open Decisions Required from PI" section still awaiting
  input (journal tier, SPURS vs RaSP, ProSST vs SaProt vs ESM-2, EVcouplings vs EVE, wet-lab
  timeline). Also not executed in the current Snakefile (still the old phase1–phase9 structure).
  Substantive, not bloat.
- **`temporary_TODO_42526.md`** (despite its name) is the live, detailed manual-retrieval task
  list — BRAKER2 SLURM commands, exact ENA/JGI download steps and gene IDs — that
  `STATUS_2026-05-12.md` and `PROGRESS.md` both directly reference and depend on. None of this
  remediation pass's five tasks resolved these manual-retrieval gaps (Tarnita 2023 genomes,
  *N. mirabilis* ENA proteins, *U. gibba* JGI proteins remain TODO in `config/enzyme_families.yaml`
  exactly as before). Actively load-bearing.
- **`PROGRESS.md`** predates and is explicitly self-superseded by `STATUS_2026-05-12.md` (per
  that file's own header: "Supersedes: PROGRESS.md"), but still contains operational detail
  (exact MAFFT/trimAl/IQ-TREE commands, Phase 2 test gates) not fully duplicated in the newer
  file. Left as-is rather than deleted, since removing it would lose that operational detail
  without anywhere else to fold it.
- **`FOLDX_REVIEW_2026-05-12.md`** and **`FOLDX_ALTERNATIVES_AF3_2026-05-12.md`** are a real,
  detailed literature review (benchmark tables against Megascale/S669, specific method
  recommendations: SPURS, CatOpt) backing the same still-open ESM-2/FoldX-stack pivot decision
  above. Not bloat.

**Recommendation for the PI, not acted on here:** these six files collectively describe one
coherent, apparently-abandoned-or-paused architectural pivot from 2026-05-12. Worth an explicit
PI decision — resume it, formally shelve it with a note, or fold its still-relevant parts into
`TASKS.md`/`CLAUDE.md` and archive the rest — rather than leaving it in an ambiguous state where
it's unclear if it's active or historical. This is a separate, larger decision from anything in
this remediation pass and was not made unilaterally here.

**Update 2026-08-24, per PI request:** re-read `TASKS.md` and found this was a more accurate
statement than intended above — `TASKS.md` already substantially tracks the pivot as Tier 0/1
checklist items (T0.1 ProSST/SaProt swap, T0.2 new genomes, T0.3 FEP cap, T1.3b DL ΔΔG ensemble),
just with every checkbox still unchecked. What was genuinely missing: the Snakefile Section A–G
reorganization itself (`PROJECT_RESTRUCTURE_2026-05-12.md` §2/§5/§6) had no corresponding TASKS.md
entry, and the 5 "Open Decisions Required from PI" (§8) were nowhere surfaced as an explicit gate.
Added both: a new `T0.0` task for the Snakefile reorganization, and a blocking-decisions checklist
at the top of Tier 0, plus an explicit staleness flag (this tier has sat open, unexecuted, for
over three months as of this remediation pass). See `TASKS.md` Tier 0 header.

## 6. Final consistency sweep (repo-wide, excluding `audit/`, `.git/`, `results/`)

- `A01.073`: zero live (uncorrected) hits.
- `G03.001`: zero live hits outside one harmless URL-example comment in
  `config/enzyme_families.yaml`.
- `Tiew` / `Tiew & Goh`: zero live hits after the fixes in §3–4 above — every remaining mention
  repo-wide is inside an explanatory "this was corrected" comment.
- Bare `chitinases_gh19` (as a standalone family key, not `_class_iv`/`_class_i`): zero hits.
- `nepenthesins`-as-whole-A1B-group: not separately re-swept in this pass (already covered by
  T3/T4's own grep checks); no new hits surfaced incidentally during this task's reading.

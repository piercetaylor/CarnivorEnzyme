# Audit Changelog

## 2026-08-22 — T1: Fix detect_convergence.py species-metadata loader

Fixed critical bug in `workflow/scripts/detect_convergence.py` where `_load_species_meta()` was reading from a non-existent `species:` key, always returning zero species and silently breaking convergence detection. Loader now correctly iterates both `carnivorous_species:` and `outgroup_species:` sections, infers `carnivorous` status from section membership, and produces a dict of 23 species (18 carnivorous + 5 outgroups).

See audit/01_detect_convergence_loader_fix.md

## 2026-08-22 — T2: Reassign Caryophyllales carnivory origins

GATE FAILED — see audit doc. Retrieved and quoted Fleck & Jobson (2023, *Plants* 12:3356) directly: it states carnivory within Caryophyllales "may have evolved once," treating Nepenthaceae and Droseraceae as sister clades from a single common ancestor rather than two independent origins — a different count than the hypothesized "2 origins" this task was gated on, and Fleischmann et al. (2018) could not be retrieved/quoted from primary text. `config/species.yaml` was left unchanged.

See audit/02_carnivory_origin_reassignment.md

## 2026-08-22 — T2 follow-up: Fleischmann 2018 re-check + edit applied

Second retrieval attempt on Fleischmann et al. (2018) still could not obtain a single clean primary-source fetch (PDF mirrors failed via TLS mismatch and corrupted/image-scan text; Google Books exposes no chapter text), but five independent WebSearch queries converged on consistent, specific technical content ("once (or twice) in the non-core Caryophyllales (i.e., Nepenthales)," with losses in Ancistrocladaceae and Dioncophyllum/Habropetalum) that is commensurable with and agrees with Fleck & Jobson (2023) — same five-family clade, same single-origin-plus-losses framing, just under an alternate ordinal name ("Nepenthales" vs. APG "Caryophyllales"). Not a genuine conflict. Per PI direction, edited `config/species.yaml`: collapsed all 9 Caryophyllales carnivorous species (Nepenthes×3, Drosera×4, Dionaea, Aldrovanda) to a single `carnivory_origin: 1`; renumbered Oxalidales 4→2, Sarraceniaceae 5→3, Lentibulariaceae 6→4, Byblidaceae 7→5; rewrote the header comment block to remove the prior internal contradiction and cite both sources with quotes; added inline comments at Dionaea/Aldrovanda explaining the shared origin.

See audit/02_carnivory_origin_reassignment.md

## 2026-08-22 — T3: MEROPS correction, nepenthesins/A1B restructure, neprosin re-scope

GATES PASSED (both Task 1 and Task 2) — verified directly against live MEROPS pages
(ebi.ac.uk/merops), independently cross-checked via WebSearch, plus Butts et al. 2016/2020 and
a 2018 PeerJ digestive-enzyme review. Findings: `merops: A01.073` on the old combined
`nepenthesins` entry was wrong (A01.073 = rice nucellin, *Oryza sativa*, not nepenthesin); the
correct nepenthesin holotype is `A01.040` (*Nepenthes gracilis*, UniProt Q766C3) — confirmed as
the ONLY carnivorous-plant holotype in MEROPS family A1. `merops: G03.001` on `neprosins` was
wrong (G03.001 = strawberry mottle virus glutamic peptidase, unrelated); neprosin's actual
MEROPS-official code is `U74.001` (family U74, "unknown catalytic type," holotype "neprosin
(*Nepenthes ventrata*)," established 2017-01-06) — Tiew & Goh 2022's proposed glutamic-peptidase
"G3" reclassification is not yet MEROPS-adopted (this specific sub-claim rests on
WebSearch-paraphrased snippets, not a verbatim primary quote, since ScienceDirect/bioRxiv fetches
were blocked/rate-limited; the U74.001 MEROPS-side fact itself is solidly, directly confirmed).
Droserasin/dionain-AP/Cephalotus-AP/Sarracenia-AP have no MEROPS holotype of their own —
classified only "A1B-like by homology" (confirmed via direct MEROPS family-A1 fetch).

Applied: split `tier1.nepenthesins` into Nepenthes-only `nepenthesins` (merops A01.040) and new
`tier1.aspartic_proteases_a1b_homology` (merops null, homology-only A1B group: droserasin,
dionain-AP, Cephalotus AP, Sarracenia AP + duplicated outgroups); fixed `neprosins` merops to
U74.001; moved `neprosins` from `tier1` into a new top-level `methods_benchmark` section
(confirmed via direct code read that `Snakefile`'s `FAMILIES = list(data["tier1"].keys())`
correctly drops neprosins and picks up the new family, and that no rule/script hardcodes
`FAMILIES`/`tier1` independently). All original accessions, TODOs, and comments preserved — no
data deleted. Annotated `non_carnivorous_homologs.a1b_aspartic_protease_outgroups` entries with
`supports_family` cross-references to both split entries.

Flagged, not resolved: splitting `nepenthesins` into two tier1 keys means two separate
alignments/trees, so reproducing Fukushima et al. 2017's actual cross-lineage A1B convergence
analysis (which included Cephalotus AP alongside Nepenthes) needs a downstream merge step not
covered by this config; `config/config.yaml`'s `fep.substrates_for_fep` still has a stale
`neprosins` key and lacks one for `aspartic_proteases_a1b_homology` (harmless — `run_fep.py`
not yet built) — left untouched as out of this task's explicit file scope.

See audit/03_merops_restructure_and_neprosin_rescope.md

## 2026-08-22 — T4: Split GH19 chitinases into Class I vs Class IV

GATE PASSED. Verified directly against Yoneda et al. 2025 (*FEBS Open Bio* 15:1930-1944, PDB
9JTR) that Class I GH19 chitinases are genuinely defined by an intact N-terminal chitin-binding
(hevein/CBM18) domain, and that D. adelae's BAW35424.1 is explicitly stated to be Class I
("belongs to class I of the GH19 family"). Independently confirmed D. binata's QHN63861.1
(DbChitI-3) as Class I via Planta 2025 (PMC11725546: "extracellular class I chitinase"). The
GH19 Engineering Database (PMC8547705) confirms Class IV instead carries a CBM18 "with a
deletion" plus a "loopless" catalytic domain — a real, independently-verifiable architectural
difference from Class I, not a labeling nitpick, and one that would confound alignment-column
homology and ancestral-state reconstruction if mixed with Class IV in one alignment/tree.

Applied: split `tier1.chitinases_gh19` into `chitinases_gh19_class_iv` (the Fukushima et al.
2017 Fig. 3a cross-lineage convergence target; everything from the old entry except the two
Class I accessions) and new `tier1.chitinases_gh19_class_i` (Drosera_adelae BAW35424.1 +
Drosera_binata QHN63861.1, labeled "structural comparison only — not used for cross-lineage
convergence detection"). Replaced the Drosera_capensis "use BDB33378.1 (class I) as proxy if
class IV unavailable" framing with an honest `TODO` (DCAP_0533 exists only in an unsubmitted
draft genome, no public Class IV accession available) — matching the Aldrovanda "REMOVED"
honesty convention already used elsewhere in the file. BDB33378.1 itself moved to
`chitinases_gh19_class_i` since it is genuinely Class I. All other accessions, TODOs, and
comments preserved verbatim — no data deleted.

Confirmed via direct code read that `Snakefile`'s `FAMILIES = list(data["tier1"].keys())` picks
up both new keys automatically; no Snakefile change needed. Grepped `workflow/rules/*.smk` and
`workflow/scripts/*.py` for the literal string `chitinases_gh19` — no functional hardcoded
references found (one harmless CLI help-string example in `detect_convergence.py`). Found and
fixed two dangling references outside that grep scope that would otherwise have silently broken
after the split: `config/config.yaml` (`fep.substrates_for_fep`) and `config/substrates.yaml`
(`GlcNAc4.for_family`), both renamed `chitinases_gh19` → `chitinases_gh19_class_iv`. Left
`workflow/scripts/braker2/04_run_hmmer_scan.sh` unfixed (out of scope — a standalone manual
helper not invoked by any Snakemake rule, already independently stale re: the prior A1B split).

See audit/04_gh19_class_split.md

## 2026-08-22 — T2 checker follow-up: latent Windows encoding bug fixed

T2's opus checker flagged that `workflow/scripts/detect_convergence.py`'s `_load_species_meta()`
opened `species.yaml` with `open(species_yaml)` (no explicit encoding), which throws
`UnicodeDecodeError` on Windows (cp1252 default) given the file's non-ASCII characters (14 such
lines pre-existing, 19 after the T2 header-comment rewrite). Harmless on the Linux/HPC deployment
target but a real portability bug. Fixed directly (one-line change, `encoding="utf-8"` added) by
the orchestrating session rather than spinning up a new fixer/checker pair, since the root cause
and fix were already fully diagnosed and evidenced by the T2 checker's own reproduction of the
crash — no new verification needed beyond what's recorded in
audit/02_carnivory_origin_reassignment.md's checker section.

See audit/01_detect_convergence_loader_fix.md, audit/02_carnivory_origin_reassignment.md

## 2026-08-22 — T3 corrections applied directly (fixer agent session lost)

The fixer agent tasked with correcting T3's DISPUTED citation errors and moving `nepenthesins`
to `methods_benchmark` was interrupted mid-task (background process teardown, not a work
failure) and did not complete either follow-up. Verified via direct file inspection that neither
change had landed. Applied both directly, using the opus checker's already-independently-verified
findings from `audit/03_merops_restructure_and_neprosin_rescope.md` §V.7 as the source of truth
(no new literature re-fetch needed — those sources/quotes were already primary-verified):

1. **Neprosin citation corrected**: "Tiew & Goh 2022" → **Ting TY, Baharin A, Ramzi AB, Ng CL,
   Goh HH (2022), *Plant Physiology and Biochemistry* 183:23-35**. Reframed the claim per the
   checker's finding: the published paper says neprosin "founds a new glutamic peptidase family"
   and does NOT itself claim MEROPS "G3" membership (that wording existed only in a superseded
   bioRxiv v1 preprint). Added Oda & Wlodawer 2023 (*Biochemistry* 62:672-694) as a separate,
   narrower citation for neprosin's structural (not formal taxonomic) similarity to family G3.
2. **"Butts et al. 2020" corrected** → **Sprague-Piercy et al. 2020, *Biomolecules* 10:1069**
   (DOI 10.3390/biom10071069) for the droserasin PSI-length citation.
3. **`nepenthesins` moved** from `tier1:` to `methods_benchmark:` (alongside `neprosins`), per
   PI decision: Nepenthes-only family is now entirely within one `carnivory_origin` (post-T2
   collapse), so it cannot produce a cross-lineage convergence result under `min_lineages: 2`.
   All accessions/comments preserved verbatim in the move. Verified `aspartic_proteases_a1b_homology`
   (the other split-off entry) still legitimately spans 3 distinct carnivory_origin values under
   the current `config/species.yaml` numbering (Drosera+Dionaea=1, Cephalotus=2, Sarracenia=3) —
   remains a valid standalone tier1 convergence target. Added `experimental_pdb: null` to it.
4. **`workflow/scripts/braker2/04_run_hmmer_scan.sh` fixed**: renamed hardcoded
   `"nepenthesins"`→`"aspartic_proteases_a1b_homology"` and `"chitinases_gh19"`→
   `"chitinases_gh19_class_iv"` (the latter was also flagged, separately, by T4's audit as
   out-of-scope-but-stale; fixed here since I was already correcting the same script).
5. **Environment correction**: earlier audit docs in this session state "no working Python 3"
   in this workspace, based on `python` resolving to an MGLTools 2.7 install. Direct check found
   `python3` (a separate PATH entry) is Python 3.13 with `pyyaml`, `click`, `biopython` 1.87,
   `scipy`, `numpy`, and `requests` all installed and working — confirmed by successfully
   `yaml.safe_load`-parsing the live `config/enzyme_families.yaml` after all edits above (tier1
   keys: `chitinases_gh19_class_iv`, `chitinases_gh19_class_i`, `purple_acid_phosphatase`,
   `rnase_t2`, `aspartic_proteases_a1b_homology`; methods_benchmark keys: `nepenthesins`,
   `neprosins`; confirmed `nepenthesins` no longer in `tier1`). Missing: `ete3` (blocks running
   `detect_convergence.py` itself, which needs it for tree parsing) and `pytest`. This does NOT
   fully overturn the earlier "can't execute the pipeline here" conclusion, but it's less broken
   than previously stated — flagging for the PI in case it changes the Stage C (re-fetch) approach.

See audit/03_merops_restructure_and_neprosin_rescope.md

## 2026-08-22 — T3 final: last stray citation fixed, checker CONFIRMED

Independent checker (sonnet) verified the direct corrections above and found one remaining
miss: `config/enzyme_families.yaml`'s `methods_benchmark.neprosins.notes:` field (a sibling
field to the already-corrected `merops_note:`/key-literature block) still said "Tiew & Goh
2022." Fixed directly (one-word author correction to "Ting et al. 2022"). T3 is now fully
CONFIRMED: MEROPS codes verified against live MEROPS pages, subfamily restructuring mechanically
complete (0 accessions lost), neprosin citation corrected to Ting et al. 2022 (*Plant Physiol.
Biochem.* 183:23-35) with accurate framing (paper claims "a new" family, not formal G3
membership; Oda & Wlodawer 2023 separately notes structural-only G3 similarity), both
Nepenthes-only `nepenthesins` and `neprosins` moved to `methods_benchmark` (single-origin
problem), `aspartic_proteases_a1b_homology` confirmed to still validly span 3 carnivory origins
as a standalone tier1 target, and the `braker2/04_run_hmmer_scan.sh` helper script updated to
match. `CLAUDE.md` §10 still has the stale "Tiew & Goh 2022... 183:23-38" citation — deferred to
the T6 documentation-sync task, not fixed here.

See audit/03_merops_restructure_and_neprosin_rescope.md

## 2026-08-22 — Stage A+B coherence gate (orchestrator, not delegated)

Cross-checked all four landed tasks against each other directly against the current
`config/species.yaml`: T1's loader reads `carnivory_origin` generically (no hardcoded values,
so it correctly consumes whatever T2 produced). `aspartic_proteases_a1b_homology` (T3) spans
origins {1 (Drosera+Dionaea), 2 (Cephalotus), 3 (Sarracenia)}. `chitinases_gh19_class_iv` (T4,
landed before T2's origin collapse was finalized — re-checked here against the final numbering)
spans origins {1 (Nepenthes, Drosera_rotundifolia, Dionaea), 2 (Cephalotus), 3 (Sarracenia)}.
Both remain valid `min_lineages: 2` cross-lineage convergence targets under the final 1-origin
Caryophyllales scheme. `chitinases_gh19_class_i` (Drosera_adelae, Drosera_binata) was always
structural-comparison-only, so origin diversity doesn't apply to it. No inconsistencies found.
Stage A+B gate: PASS. Proceeding to Stage C.

## 2026-08-22 — Correction: .env was never actually empty

Earlier in this session (and baked into this plan / several audit docs) it was claimed
`.env` has empty `NCBI_EMAIL`/`NCBI_API_KEY`, based on running `grep -o "^[A-Z_]*=" .env`.
That command's `-o` flag only prints the matched pattern itself (the key name plus `=`),
not the rest of the line — so it displayed `NCBI_EMAIL=` regardless of whether a value
followed. This was a diagnostic error by the orchestrating session, not an actual gap in
the project: `.env` already has a real `NCBI_EMAIL` and `NCBI_API_KEY` configured. Using
the existing credentials for Stage C rather than substituting a different email address.

## 2026-08-24 — T6: Documentation sync + bloat trim

Updated `CLAUDE.md` (Section 3 family table, Section 4 origin correction, Section 10 citation),
`README.md`, and `TASKS.md` to match the corrected `config/enzyme_families.yaml` and
`config/species.yaml`. Fixed the same stale "Tiew & Goh 2022" neprosin misattribution (should be
Ting TY, Baharin A, Ramzi AB, Ng CL, Goh HH 2022, *Plant Physiol. Biochem.* 183:23-35) in
`literature_review/LITERATURE_FINDINGS.md`, `literature_review/MASTER_LITERATURE_REGISTRY.md`
(which also had a wrong DOI), and `ORTHOLOGY_AND_FAMILY_SCOPE_2026-05-12.md` — two of these three
files were outside the original task checklist and only surfaced via a full-repo grep sweep.
Investigated 6 root-level files flagged as likely-superseded "bloat" by an earlier audit pass
(`STATUS_2026-05-12.md`, `PROGRESS.md`, `PROJECT_RESTRUCTURE_2026-05-12.md`,
`temporary_TODO_42526.md`, `FOLDX_REVIEW_2026-05-12.md`, `FOLDX_ALTERNATIVES_AF3_2026-05-12.md`)
by reading all six in full — **none were bloat**; all are substantive, still-load-bearing
documents (a real, apparently-unexecuted ESM-2→ProSST architecture pivot proposal with an open
PI-decision section, and the live manual-retrieval TODO list still referenced by
`config/enzyme_families.yaml`'s TODO entries). None deleted; flagged for a PI decision instead.
Final repo-wide grep sweep confirms zero live (uncorrected) references to A01.073, G03.001, or
"Tiew & Goh" remain outside `audit/` and explanatory correction comments.

See audit/06_doc_sync_and_bloat_trim.md

## 2026-08-24 — T5: Live re-fetch via the real fetch_sequences.py (not WebFetch)

Given the working `python3` environment discovered above, ran the actual production
`workflow/scripts/fetch_sequences.py` against live NCBI Entrez / UniProt REST, using the
existing `.env` credentials, rather than hand-writing FASTA files via WebFetch as originally
planned — more faithful to the real pipeline (same code path, same logging, same length-filter
logic that will run on Hellbender).

**Script gap found and fixed first**: `fetch_sequences.py` only ever iterated
`config["tier1"]`, with no awareness of the new `methods_benchmark:` section T3 created —
meaning `nepenthesins` and `neprosins` would have been silently skipped entirely, and the old
stale `results/sequences/nepenthesins/` directory would have stayed orphaned. Fixed by building
the iteration dict as `{**config.get("methods_benchmark", {}), **config.get("tier1", {})}`
(workflow/scripts/fetch_sequences.py, ~line 169-181) — no key collisions between the two
sections, tier1 wins if that ever changes.

**Execution**: purged the stale `results/sequences/nepenthesins/` directory first (gitignored,
disposable pipeline output — confirmed via `git check-ignore`), then ran a full fetch across all
7 families (tier1 + methods_benchmark). Result: `fetched=92 skipped/filtered=30 failed=0`, exit
code 0.

**Contamination verified purged**: grepped the entire `results/sequences/` tree for the 4
previously-flagged bad accessions (GAV28673.1 yeast contaminant, GAV69561.1
UDP-glucuronosyltransferase, KAL9250314.1 non-digestive guard-cell isoform, XP_002267025.1
wrong-enzyme phosphatidylinositol kinase) — zero hits. `Triticum_dicoccoides.fa` now correctly
labeled (was mislabeled `Triticum_aestivum.fa` in the stale data). New `nepenthesins/` directory
has exactly the 14 sequences the corrected config specifies (3 Nepenthes species x2 + 5 outgroup
species — Arabidopsis, Solanum, Hordeum x2, Triticum_dicoccoides, Vitis).

## 2026-08-24 — T5 checker: independent verification, CONFIRMED

Independent checker re-verified T5 from scratch (no trust of the changelog entry above). Grepped
all of `results/sequences/` for the 4 flagged contaminants — zero hits. Independently re-fetched
`GAV73214.1`, `GAV73783.1`, `XP_019076990.1`, `XP_037427783.1` directly from live NCBI
`efetch.fcgi` and diffed residue-for-residue against the repo files — exact matches, correct
organisms (`Triticum dicoccoides` confirmed, no stray `Triticum_aestivum.fa` anywhere). Checked
every FASTA header in all 61 files across all 8 family directories against the
`>{Species}|{acc} family=... species=... acc=... len=...` convention — zero violations. Counted
non-TODO config accessions vs. fetched sequences for 7 families: 5 matched exactly; the 2 that
didn't (`aspartic_proteases_a1b_homology` 16→14, `purple_acid_phosphatase` 16→15) were both
confirmed to be correct, intentional drops by the `min_length_fraction=0.8` median-length filter
(reproduced the filter arithmetic by hand; independently re-fetched the dropped
`KAL6998072.1` from NCBI to confirm its exact length of 428 aa, just under its family's 431.2 aa
cutoff). Confirmed the `{**methods_benchmark, **tier1}` merge in `fetch_sequences.py` has no key
collisions and both `nepenthesins/`/`neprosins/` are fresh and non-empty. Confirmed `results/` is
fully gitignored and untracked.

One non-blocking gap found and flagged (not a contamination issue): a stale orphaned
`results/sequences/chitinases_gh19/` directory (pre-dating the T4 class split, Apr 14 mtime) was
never purged by the re-fetch — harmless (unread by any Snakemake rule, no bad data), but should be
deleted in a future pass.

See audit/05_live_refetch_verification.md

## 2026-08-24 — T9: Fix silent-failure patterns in 04_run_hmmer_scan.sh

Fixed four issues in `workflow/scripts/braker2/04_run_hmmer_scan.sh`, each verified by running
the actual bash/python logic in isolation with synthetic stand-in files, then confirmed end-to-end
against the real edited script with stubbed `hmmsearch`/`hmmfetch` binaries (no real hmmer binary
or proteome data available in this environment): (1) the zero-hits guard was unreachable under
`set -euo pipefail` because `grep -v '^#'` exits 1 on an all-comment `.tblout`, killing the script
silently before the intended "0 hits" message — added `|| true`; (2) if all 5 target proteomes are
missing the script exited 0 with a misleading success summary — added an `n_scanned` counter gated
before the `=== Summary ===` block, exits 1 with an explicit error if nothing was scanned; (3)
hmmsearch hit IDs absent from the proteome FASTA (e.g. stale `.tblout`) were silently dropped and
the printed count wrongly blamed the length filter — now reports resolved-vs-total counts and
warns explicitly on any unresolvable IDs; (4) a leftover pre-split family-name path
(`..._chitinases_gh19.faa`) in the final "Next steps" example was renamed to
`..._chitinases_gh19_class_iv.faa` (confirmed via full-file grep: zero remaining stale references).

See audit/09_hmmer_scan_script_fixes.md

## 2026-08-24 — T10: Fix config structural gaps, add consistency test

Moved `chitinases_gh19_class_i` from `tier1:` to `methods_benchmark:` in
`config/enzyme_families.yaml` (verbatim, all comments/accessions preserved) — its 3 species
(Drosera_capensis, Drosera_adelae, Drosera_binata) all share `carnivory_origin: 1`, so under
`convergence.min_lineages: 2` it could never produce a cross-lineage convergence result, and
because `Snakefile:44-49` derives `FAMILIES` from `tier1.keys()` alone, it was silently baked
into every convergence-related Snakemake target (phase2/phase3a/phase4a/phase5d/rule all) as a
guaranteed header-only empty output — the exact "zero rows, exit 0" false-positive signature
this remediation effort targets. Same treatment already applied to nepenthesins/neprosins on
2026-08-22 for the identical single-origin reason.

Found and fixed 7 species-config gaps via a real PyYAML cross-reference script (before: 7
missing; after: 0): added `Drosera_rotundifolia`, `Nepenthes_ampullaria`,
`Nepenthes_bicalcarata`, `Nepenthes_rafflesiana`, `Nepenthes_ventrata` to
`species.yaml.carnivorous_species` and `Triticum_dicoccoides` to `species.yaml.outgroup_species`
(all taxids/genome accessions WebSearch/NCBI-eutils-verified, none invented; two TODO-only
Nepenthes species given lighter-weight tier:2 placeholders per judgment call); renamed the
malformed `Nepenthes_mirabilis_PAP` accession key in `tier1.purple_acid_phosphatase` to the
correct `Nepenthes_mirabilis` (confirmed no collision — plain rename, not a merge).

Added `tests/test_config_consistency.py` (5 test functions covering the requested (a)-(d) plus
one extra direct regression guard) and manually ran every assertion's logic via `python3 -c`
against the real config files (pytest itself isn't installed in this environment) — all 5 pass,
including a round-trip check of `detect_convergence.py`'s `_seq_id_to_species()` against every
species.yaml name, which also re-confirms the species-name-vs-code keying fix from
`audit/07_detect_convergence_deep_fix.md` is still correctly in place.

See audit/10_config_gaps_and_consistency_test.md

## 2026-08-24 — T8: Harden fetch_sequences.py error handling

Fixed 5 silent-failure modes in `workflow/scripts/fetch_sequences.py` found in code review: (1)
malformed non-list `accessions:` entries now log ERROR + increment `n_failed` instead of vanishing
silently; (2) `--families` now validates requested keys against the real family set and raises
`click.BadParameter` on a stale/typo'd key instead of silently running an empty fetch; (3) fetched
records are now checked for an accession mismatch *before* being relabeled with the requested
accession string, so a wrong-protein NCBI/UniProt resolution is caught and excluded rather than
silently mislabeled as evidence-destroying (self-corrected a bug in the initial mismatch-check
implementation during verification: UniProt/Swiss-Prot `sp|ACCESSION|NAME` records need field
index 1, not the last field, confirmed via live fetch of `P42210`); (4) `tier1`/`methods_benchmark`
key collisions now raise `ValueError` instead of silently resolving via dict-merge precedence; (5)
families missing `expected_length_aa` now log a WARNING that length QC is disabled (currently
prophylactic — all 7 fetched families have the key). All 5 fixes verified by real CLI execution
against the real config and/or scratch throwaway YAMLs/output dirs. See
audit/08_fetch_sequences_hardening.md.

## 2026-08-24 — T7: Deep fix for detect_convergence.py (species keying, .state parsing, fabricated fallbacks)

Fixed 4 independent defects in `workflow/scripts/detect_convergence.py`, EACH of which alone caused
the script to emit zero convergent sites with no error (the prior remediation `96ab469` fixed none
of them — verified by re-running the pre-fix code). (1) `_load_species_meta()` keyed its dict by the
4-letter `code` field while every tree leaf/FASTA header is `{Species_name}|{accession}`, so
`_seq_id_to_species()` matched nothing — now registers each species under BOTH the species-name and
code keys; also subsumes HIGH-3 by raising on a carnivorous species missing `carnivory_origin`
instead of silently defaulting it to `0` (which the detector then drops). Verified against the real
`config/species.yaml`: all 18 carnivorous entries have an origin, all 5 origins present, no
ValueError — no latent data problem. (2) `_parse_state_file()` read `parts[-1]` (a `p_V` float
string) as the ancestral amino acid and never used the header it parsed; its `AA_COLS.index(aa) + 2`
offset was ALSO off by one (real header is `Node Site State p_A …`, so `p_A` is index 3, not 2), and
a bare `except` swallowed the resulting ValueError into `posterior = 0.0`. Now resolves `State` and
`p_<State>` by header name, lets a missing column propagate, and raises if 0 internal nodes parse.
Format verified against the IQ-TREE Command Reference; confirmed IQ-TREE's AA order
`ARNDCQEGHILKMFPSTWYV` does match the file's `AA_COLS` constant. (3) `_load_tree()` treated
IQ-TREE `-bb/-alrt` support labels (`85.7/97`) as node names, so they never matched `.state` NodeN
labels — support labels are now detected and moved aside, plus a defensive tree-vs-.state node-name
overlap check in `main()` that aborts loudly instead of yielding zero rows. (4)
`_get_carnivorous_branch_lengths()` returned a fabricated `1.0` and `_background_substitution_rate()`
a magic `0.01` when they had no data; both feed `mu_per_site` for the Poisson test, so the pipeline
was producing publication-shaped p-values from invented inputs. Both now raise. All fixes verified by
real `python3` execution against the real config and a synthetic IQ-TREE-format `.state` file;
end-to-end run proves old code = 0 rows and fixed code = 1 row on identical input (`ete3` is not
installed, so `_load_tree` itself is trace-verified only — its extracted `_is_support_label`
predicate and the `main()` overlap check were executed).

See audit/07_detect_convergence_deep_fix.md

## 2026-08-24 — T13: Distinguish real I/O errors from legitimate zero-hits in 04_run_hmmer_scan.sh

Replaced the blanket `|| true` around the `grep -v '^#' "$hit_tbl" | awk '{print $1}' | sort -u`
pipeline in `workflow/scripts/braker2/04_run_hmmer_scan.sh` with logic that isolates the `grep`
call and inspects its exit code: `grep_rc == 1` (no non-comment lines — a legitimate empty hmmer
result) is still tolerated and reported as "0 hits", but `grep_rc >= 2` (missing/unreadable file,
real I/O error) now aborts the script with a descriptive `ERROR:` message instead of being
silently misreported as "0 hits". `awk`/`sort` are no longer swallowed by the same blanket guard.
Verified with 3 scratch tests against synthetic `.tblout` files (comment-only file, nonexistent
file, real hit lines) — all three transcripts in audit/13_hmmer_scan_error_distinction.md.

See audit/13_hmmer_scan_error_distinction.md

## 2026-08-24 — T12: fetch_sequences.py — make silently-empty species/families visible

Fixed two remaining silent-failure patterns in `workflow/scripts/fetch_sequences.py` left over
from T8's hardening pass: (1) a species whose every fetched record was excluded by the per-family
length filter wrote no FASTA and left no trace beyond a generic `n_skipped` increment that never
affects the exit code — now logs a `ZERO SURVIVING` WARNING per occurrence and is tallied in a new
`n_species_all_filtered` counter reported in the final summary; (2) a family whose every accession
is still `TODO` (or, extended beyond the original ask, whose every species gets entirely
length-filtered) logged a single WARNING and moved on — now logs `ZERO OUTPUT` at ERROR level and
the family key is collected into a `families_zero_output` list, unconditionally reported at ERROR
level in the final summary. Neither condition forces `sys.exit(1)` by itself (only real fetch
failures do): both can be legitimate, already-documented outcomes (e.g. `enzyme_families.yaml`'s
own comment on Dionain-3, `A0A0E3M338`, flags it as "very likely a partial/truncated gene model...
Flag for length filter in pipeline" — forcing a hard failure on exactly the outcome the filter is
designed to produce would be wrong), so visibility (not exit-code severity) is the fix.

Checked the real config first: none of the 7 currently-fetched families (`tier1` +
`methods_benchmark`) has zero non-TODO accessions — Fix 2's family-level TODO case is prophylactic
today. But a live run against the real config with real NCBI/UniProt fetches found **2 real,
previously-invisible occurrences of Fix 1's condition**: `purple_acid_phosphatase/
Sarracenia_purpurea` (`KAL6998072.1`, 428 aa, just under the family's 431 aa cutoff — the same
drop the T5 checker had found only via manual arithmetic after the fact) and
`aspartic_proteases_a1b_homology/Dionaea_muscipula` (both of its two accessions, 352 aa and 300
aa, fall under the family's 365 aa cutoff). Both fixes additionally verified against a disposable
synthetic YAML built to force both conditions deterministically (`species_all_filtered=1`,
`families_zero_output=1`, both named explicitly in the summary). All commands run for real against
either the live config or scratch YAMLs/output dirs; `results/sequences/` was never touched.

See audit/12_fetch_sequences_visibility_fixes.md

## 2026-08-24 — T14: Close test blind spots, fix CLAUDE.md/README tier-list drift

Closed two test-coverage gaps in `tests/test_config_consistency.py` found by review: (a) added
`test_carnivorous_and_outgroup_species_disjoint`, guarding against a species listed in both
`carnivorous_species:` and `outgroup_species:` in `config/species.yaml` — `detect_convergence.py`'s
`_load_species_meta()` loads carnivorous first then outgroups second, so the second silently
overwrites the first, reclassifying a carnivorous species as non-carnivorous
(`carnivory_origin=0`); the pre-existing union-based `_known_species_keys()` helper cannot detect
this by construction. Verified: duplicated `Cephalotus_follicularis` into a scratch copy of
`species.yaml` and confirmed all 5 pre-existing tests still passed (the blind spot is real), the
new test fails against that scratch file, and passes against the real file (no duplication exists
today). (b) Added `test_species_codes_are_unique` for the same silent-overwrite class keyed on
`code:` instead of species name — verified the same way with a scratch duplicate `code`. (c) Added
`test_tier1_families_span_min_lineages_carnivory_origins`, asserting every `tier1:` family's
non-TODO accessions span `>= config.yaml`'s live `convergence.min_lineages` (2) distinct
`carnivory_origin` values — the exact invariant that motivated moving
`chitinases_gh19_class_i`/`nepenthesins`/`neprosins` out of `tier1:` in T3/T10. Verified against
the real config (all 4 `tier1:` families pass, `rnase_t2` right at the n=2 threshold) and against
a scratch config with `aspartic_proteases_a1b_homology` reduced to a single origin (test
correctly fails). (d) Fixed `test_seq_id_to_species_roundtrip_on_species_names`, whose docstring
claimed to guard `_load_species_meta()`'s species-name-vs-code keying fix but actually built its
`known_species` set independently from `species.yaml` rather than from the loader's real output —
verified by reverting the keying fix in a scratch copy of `detect_convergence.py` and confirming
the old test body still passed (0 failures) against the broken loader while the fixed body (now
sourcing `known_species` from `module._load_species_meta(...).keys()`) correctly failed (29
failures, one per species). All 8 tests in the file (5 original, one docstring-fixed, 3 new) pass
against the real, unmodified config files, run manually via `python3` since `pytest` is still not
installed in this environment.

Fixed a doc-drift violation: `CLAUDE.md` §3 still listed `chitinases_gh19_class_i` as Tier 1 item
#2 — directly contradicting its own "this list must match `config/enzyme_families.yaml`'s actual
`tier1:` keys exactly" rule, left stale since the T10 move to `methods_benchmark:`. Renumbered
Tier 1 to its real 4 items and moved `chitinases_gh19_class_i` into the Methods Benchmark
subsection with the same single-origin framing as its `nepenthesins`/`neprosins` siblings.
`README.md` had the matching gap (GH19 Class I row missing the `methods_benchmark.` prefix its
siblings have, plus a stale "retained tier1" lead-in sentence) — fixed the same way. Verified
programmatically (not by eyeballing): a `python3` script PyYAML-parses the real
`tier1`/`methods_benchmark` key sets and confirms both documents' family lists now match exactly.

See audit/14_test_blindspots_and_doc_drift.md

## 2026-08-24 — T11: Architectural fix for ancestral-reconstruction/tree-rooting mismatch
Split IQ-TREE into two passes so ancestral state reconstruction runs on the SAME topology the
pipeline consumes downstream, eliminating the node-correspondence problem rather than patching it.
Pass 1 (`infer_tree`) now does topology search only (`-bb`/`-alrt`, no `-ancestral`); `root_tree.py`
additionally emits a label-stripped rooted tree and the matched outgroup tip labels; a new
`ancestral_reconstruction` rule runs `run_ancestral.py`, which reuses pass 1's ModelFinder verdict
and calls `iqtree3 -m <model> -te <rooted.plain> -o <outgroups> --ancestral`, then verifies that the
emitted `.asr.treefile` and `.asr.state` describe exactly the same node set and that the topology
was unchanged. `convergence.smk` and `ancestral_structure.smk` now consume that self-consistent
pair. Added `wildcard_constraints` to the Snakefile (`{family}.asr.treefile` was otherwise also
matchable by `infer_tree` with `family="…asr"`).

In `detect_convergence.py`: deleted `_is_support_label()` and the `NodeN` counter in `_load_tree()`
— the script no longer invents any node name — and replaced the set-overlap heuristic with
node-by-node checks that raise if any internal node, or any leaf's direct parent, fails to resolve
in `.state`, or if any tree leaf is missing from the alignment. Folded in the requested visibility
fix: unmapped leaves are now counted and listed at WARNING (was one `debug` line), per-origin leaf
counts are logged at INFO, and the run raises when fewer distinct carnivory origins are represented
than `--min-lineages` (`--allow-insufficient-origins` to override for single-origin
methods_benchmark families) — origins 2 and 5 rest on a single species each, so one unparsed leaf
label could silently delete a whole independent origin.

Verified by real execution, not documentation reading: downloaded and ran the official IQ-TREE
3.1.3 Windows binary against a 9-taxon AliSim-simulated alignment using this pipeline's real
`{Species_name}|{accession}` labels. Reproduced the shipped bug end-to-end (overlap check passed on
one coincidental `Node1` match while 6/7 nodes resolved to nothing, exit 0, empty output); confirmed
`-alrt` + `--ancestral` emits compound `Node3/99.7/100` labels AND drops the root label while
`.state` still lists it; confirmed IQ-TREE logs `rooted tree` then silently unroots a rooted `-te`
input, dissolving the root vertex; confirmed `-B` is rejected with `-te` ("Ultrafast bootstrap does
not work with -fast, -te or -n option"); confirmed marginal ASR is root-invariant (max |Δp| = 1e-5
over 7 nodes × 300 sites). Then ran the fixed chain end-to-end (7/7 nodes and 9/9 leaf ancestors
resolve, vs 1/7 and 2/9 before), a positive control recovering a planted 3-origin P→W substitution
at exactly the planted position with exactly the planted origin set, four negative controls all
raising, `snakemake -n` showing correct 4-rule chaining with no ambiguity, and the existing 13-test
suite passing. NOT verified: Snakemake *execution* of the new rules (its shell layer fails on this
Windows laptop for even a trivial `echo` rule) and behaviour at real alignment scale — both need a
Hellbender run. Also surfaced a pre-existing, unrelated `ChildIOException` in `retrieve.smk` that
currently blocks any full-workflow Snakemake run; left unfixed and flagged for follow-up.

See audit/11_ancestral_reconstruction_architecture_fix.md

---

## 15 — Stability predictor audit (2026-09-01)

Audited the ΔΔG/stability axis against the literature and against first-hand execution. **Reversed
the 2026-08-25 "SPURS only, FoldX dropped entirely" decision** (PROJECT_PLAN.md §2.2): FoldX 5.1 is
reinstated as the primary stability predictor with RaSP as a single structure-sensitivity
cross-check, SPURS demoted to optional, the stability axis re-scoped from thresholded kcal/mol to
within-family percentile rank, the FoldX pH sweep demoted to a supplementary check, CpHMD downgraded
to exploratory, and Phase 6 electrostatics promoted to the project's primary pH axis. Cut CatOpt,
EpHod, the PypKa-for-PROPKA swap, and the multi-predictor "consensus ensemble" outright.

Verified by execution, not documentation reading: FoldX 5.1 (installed locally) was run end-to-end
on 7ZVA chain A — RepairPDB, then PositionScan over 87 mutations at 16 positions stratified by
Shrake–Rupley relative SASA into buried/exposed × titratable/non-titratable, at pH 2.5, 5.0 and 7.0.
Result: mean |ΔΔG| is 6.1–7.1 kcal/mol at buried positions but **0.39–0.42 at exposed positions**,
with 0/43 exposed mutations exceeding the published ±2.9 kcal/mol 95% prediction interval and 2/43
exceeding the 1.25 kcal/mol RMSE — so `config.yaml`'s 1.0 kcal/mol threshold is degenerate at
exactly the sites Fukushima 2017 reports the convergent substitutions occupy. The pH sweep fares
worse: the pH 2.5 → 5.0 shift at exposed sites averages **0.12 kcal/mol** (max 0.50), ~10× below the
tool's own RMSE, none of 43 clearing it.

Also established from the primary text that Fukushima et al. 2017 already reported the convergent
positions "do not overlap with or cluster around catalytically essential amino acids… they tend to
be located at exposed positions… likely to change their interactions with other molecules in
solution, **rather than changing protein conformation**" — i.e. the source paper describes the
observable ΔΔG-of-folding does not measure, and points instead at surface electrostatics.

Two independent citation-verification passes found **five fabricated references** in the doc set:
"Botte et al." (FoldX 2025 — real: Delgado, Reche, Cianferoni et al.), "Cai" and "Cao" (SPURS — real:
Li & Luo), "Ressa" (JanusDDG — real: Barducci et al.), and two invented GROMACS CpHMD references
sharing one bogus page number (JCTC 18:6320, actually an unrelated IDP-dynamics paper). In every
case the DOI was correct and the author invented. Also corrected: the "MAVISp March 2026 / FoldX PCC
0.40 on Megascale" claim exists in no source; CLAUDE.md's 0.711/1.238 FoldX figures describe a
median-of-five-runs protocol, not the force-field revision; GROMACS constant-pH MD is not in any
release through 2026.1 and its fork's README asks users not to publish from it.

Six superseded root-level docs moved to `archive/` with headers. NOT verified: nothing was run on
Hellbender or on AF3-predicted structures — the SASA-stratified scan used the 7ZVA crystal, and the
same check must be repeated on real AF3 output once Phase 3 lands (TODO T1.2).

See audit/15_stability_predictor_audit.md

---

## 16 — Snakemake DAG ChildIOException fix (2026-09-02)

Cleared TODO T0.2, the last local blocker in front of the Hellbender Tier 0 verification run.
`retrieve.smk` declared `fetch_all_sequences` → `directory("results/sequences/")` and
`combine_family_sequences` → `results/sequences/{family}/combined.fa`; Snakemake rejects a DAG where
one rule's output is a child of another's, so **every dry run past Phase 1 aborted before scheduling
a single job**. Fixed by moving the derived concatenation out of the checkpoint's subtree to
`results/family_fasta/{family}.combined.fa`.

Second defect found and fixed in the same pass: `esm2.smk` declared an input
`results/sequences/{family}/all_sequences.fa` that **no rule in the workflow produces** — the file
`combine_family_sequences` actually writes is `combined.fa`. Phase 3A would have failed with
MissingInputException the moment the ChildIOException stopped masking it.

Third defect found and **deliberately left unfixed**: `snakemake -n` with no target resolves to
`rule phase1` (Snakefile:83), not `rule all` (Snakefile:141), because Snakemake uses the first rule
in the file. The Snakefile's own usage docstring claims the bare invocation runs the full pipeline;
it downloads sequences and exits 0. Changing a workflow's default target is a behavioural decision,
so it was filed as TODO T0.5.5 rather than folded into a bugfix.

Verified by execution, not inspection (Snakemake 9.25.2): the T0.1 target now builds an 8-job DAG in
the correct order — fetch → combine → align → trim → infer_tree → root_tree →
ancestral_reconstruction → detect_convergence, the exact four-rule Phase-2 chain audit/11 depends on.
`phase2` 26 jobs, `phase3a` 34 jobs, `phase1` 2 jobs, all clean. `phase4a`, `phase5d` and `all` now
fail only on `results/structures/` and `results/atlas/*` from the stubbed `predict_structure.smk`
and `integrate.smk` — expected until Tier 2. 13/13 tests still pass.

NOT verified: nothing ran on Hellbender, and no job was executed at all — `-n` builds the DAG and
says nothing about whether fetch_sequences.py, MAFFT or IQ-TREE succeed on real input. Also noted
for the cluster run: the dry run schedules `fetch_all_sequences` to re-run on an mtime trigger from
`config/enzyme_families.yaml`, so an rsynced `results/sequences/` needs `snakemake --touch` first or
NCBI access on the compute node.

See audit/16_snakemake_dag_childio_fix.md

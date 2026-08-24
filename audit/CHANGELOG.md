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

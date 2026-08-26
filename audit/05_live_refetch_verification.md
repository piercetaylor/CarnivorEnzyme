# Audit 05: Independent Verification of Live Re-fetch (T5)

**Date:** 2026-08-24
**Task:** Independently verify the claims in `audit/CHANGELOG.md`'s "T5: Live re-fetch via the
real `fetch_sequences.py`" entry — that a live NCBI Entrez / UniProt REST re-fetch, run after
fixing a `tier1`-only iteration bug in `workflow/scripts/fetch_sequences.py`, produced clean,
correct, complete `results/sequences/` output. This is a from-scratch independent check: sources
were re-fetched directly, files were re-read fresh, nothing from the T5 changelog entry was
trusted without reproduction.

---

## 1. Contamination check — 4 previously-flagged bad accessions

Grepped the entire `results/sequences/` tree (all families, all species files) for each of the
four accessions flagged as contaminants/mis-assignments in earlier audit rounds:

| Accession | Previously wrongly used as | Hits found |
|---|---|---|
| `GAV28673.1` | Cephalotus AP (actually *Pichia* yeast contaminant) | **0** |
| `GAV69561.1` | Cephalotus AP (actually UDP-glucuronosyltransferase) | **0** |
| `KAL9250314.1` | *Drosera capensis* nepenthesin (actually non-digestive guard-cell isoform) | **0** |
| `XP_002267025.1` | *Vitis vinifera* outgroup (actually phosphatidylinositol 4-kinase gamma 5) | **0** |

**Zero hits for all four.** Confirmed clean.

---

## 2. Independent re-fetch of replacement accessions from live NCBI

Fetched all four directly via `efetch.fcgi` (`db=protein&rettype=fasta&retmode=text`) and diffed
against what is actually written in the repo's FASTA files (residue-for-residue, not just header
text).

| Accession | NCBI organism / annotation | Sequence match vs. repo file |
|---|---|---|
| `GAV73214.1` | *Cephalotus follicularis*, "Asp domain-containing protein/SapB_2 domain-containing protein/SapB_1 domain-containing protein" | **Exact match**, `results/sequences/aspartic_proteases_a1b_homology/Cephalotus_follicularis.fa` (517 aa) |
| `GAV73783.1` | Same annotation, *C. follicularis* | **Exact match**, same file (508 aa) |
| `XP_019076990.1` | *Vitis vinifera*, "probable aspartyl protease At4g16563" | **Exact match**, `results/sequences/nepenthesins/Vitis_vinifera.fa` (458 aa) |
| `XP_037427783.1` | ***Triticum dicoccoides***, "aspartic proteinase nepenthesin-2-like" | **Exact match**, `results/sequences/nepenthesins/Triticum_dicoccoides.fa` (469 aa). Organism confirmed **dicoccoides**, not aestivum. |

The SapB_1/SapB_2 domains on the two Cephalotus accessions are the Plant-Specific Insert (PSI)
that defines A1B aspartic proteases — annotation is consistent with the family assignment, not
just the accession number. Confirmed **no** stray `Triticum_aestivum.fa` exists anywhere under
`results/sequences/` (repo-wide filename search, zero hits); `Triticum_dicoccoides.fa` exists in
both `nepenthesins/` and `aspartic_proteases_a1b_homology/` as expected (duplicated outgroup, per
the T3 split).

---

## 3. Header format consistency

Spot-checked headers in `purple_acid_phosphatase`, `rnase_t2`, `aspartic_proteases_a1b_homology`,
`chitinases_gh19_class_iv`, and `neprosins` (five different families, not just the
originally-flagged `nepenthesins`) — all follow
`>{Species_name}|{accession} family={family_key} species={species_name} acc={accession} len={length}`.

Then ran a regex check across **every** FASTA header in the entire `results/sequences/` tree
(all 61 files) for lines *not* matching that exact pattern: **zero non-conforming headers found.**
Convention held for the whole re-fetch, not just one family.

---

## 4. Completeness against the current `config/enzyme_families.yaml`

Parsed the live config with PyYAML and counted non-TODO accessions per family (merging
`methods_benchmark` + `tier1`, matching the script's own iteration logic), then counted actual
FASTA records written per family directory:

| Family | Non-TODO accessions in config | Sequences fetched | Match? |
|---|---|---|---|
| `nepenthesins` | 14 | 14 | Yes |
| `aspartic_proteases_a1b_homology` | 16 | 14 | **2 short**, investigated below |
| `rnase_t2` | 20 | 20 | Yes |
| `chitinases_gh19_class_iv` | 18 | 18 | Yes |
| `chitinases_gh19_class_i` | 3 | 3 | Yes |
| `neprosins` | 5 | 5 | Yes |
| `purple_acid_phosphatase` | 16 | 15 | **1 short**, investigated below |

No `/tmp/fetch_run2.log` was present in this environment to consult, so both discrepancies were
independently investigated by reproducing the script's own median-length-filter arithmetic
(`workflow/scripts/fetch_sequences.py` lines 242–266: `cutoff = min_length_fraction * median(all
fetched lengths in family)`, default `min_length_fraction=0.8`):

- **`aspartic_proteases_a1b_homology`** — missing `A0A0E3GLN3` (Dionain-1, 352 aa) and
  `A0A0E3M338` (Dionain-3, 300 aa), both *Dionaea muscipula*. All 16 fetched lengths (14 present +
  these 2): `[300, 352, 393, 432, 440, 442, 449, 455, 458, 461, 469, 485, 486, 508, 508, 517]`.
  Median of 16 = 456.5; cutoff = 0.8 × 456.5 = **365.2 aa**. Both 300 and 352 fall below cutoff →
  correctly filtered. This matches the config's own inline comments on these two accessions
  ("352 aa is shorter than canonical nepenthesins... may be partial prediction"; "300 aa is at the
  lower edge... Very likely a partial/truncated gene model. Flag for length filter in pipeline") —
  the config author anticipated exactly this outcome.
- **`purple_acid_phosphatase`** — missing `KAL6998072.1` (*Sarracenia purpurea*, "putative purple
  acid phosphatase 20"). Independently re-fetched from NCBI: confirmed exact length **428 aa**.
  All 16 fetched lengths: `[428, 461, 463, 464, 469, 470, 472, 532, 546, 611, 618, 638, 640, 663,
  664, 679]`. Median = 539; cutoff = 0.8 × 539 = **431.2 aa**. 428 < 431.2 → correctly filtered
  (a narrow, legitimate miss — 3.2 aa under threshold).

Both discrepancies are the length filter working as designed on genuine short-outlier sequences,
not a fetch or contamination problem. Not a defect.

---

## 5. `methods_benchmark` merge fix — mechanical verification

Read `workflow/scripts/fetch_sequences.py` lines 163–199 directly. Confirmed:

```python
all_families: dict = {**config.get("methods_benchmark", {}), **config.get("tier1", {})}
```

`methods_benchmark` keys (`nepenthesins`, `neprosins`) and `tier1` keys (`chitinases_gh19_class_iv`,
`chitinases_gh19_class_i`, `purple_acid_phosphatase`, `rnase_t2`, `aspartic_proteases_a1b_homology`)
are disjoint sets — no collisions, confirmed by direct enumeration of both sections in the live
YAML — so the "`tier1` wins on collision" comment is currently inert but harmless. The loop then
iterates `all_families.items()` uniformly for both sections, `family_dir = output_root /
family_key`.

Confirmed both `results/sequences/nepenthesins/` (8 files, 14 sequences, all with fresh Aug 24
10:41 mtimes) and `results/sequences/neprosins/` (3 files, 5 sequences, same fresh mtimes) exist,
are non-empty, and contain valid FASTA content with correctly formatted headers. Fix works as
claimed.

---

## 6. Disposability of `results/` — git status

```
.gitignore:2:results/     results/sequences/
.gitignore:2:results/     results/
```

`git check-ignore -v` confirms `results/sequences/` and `results/` are ignored via `.gitignore`
line 2. `git status --short` shows nothing under `results/` as tracked or untracked (only
`config/*.yaml`, `README.md`, `workflow/scripts/*.py`, `workflow/scripts/braker2/*.sh`, and the
new `audit/` directory appear). `git ls-files results/` returns nothing — no file under `results/`
has ever been committed. Purging and regenerating this directory was appropriately treated as
disposable pipeline output, requiring no extra git caution.

---

## Additional finding (not one of the 6 checklist items, flagged for completeness)

**Stale orphaned directory not cleaned up by the re-fetch:** `results/sequences/chitinases_gh19/`
(8 files, dated **Apr 14 17:58** — pre-dates both the T4 class-split and the Aug 24 re-fetch)
still exists on disk. It corresponds to the *old*, pre-split `chitinases_gh19` config key, which
no longer exists anywhere in the current `config/enzyme_families.yaml` (replaced by
`chitinases_gh19_class_iv` / `chitinases_gh19_class_i` in T4). The T5 changelog entry only
describes purging the stale `nepenthesins/` directory before the run; it does not mention this
one, and indeed `fetch_sequences.py` has no cleanup step — it only ever *writes* to
`family_dir = output_root / family_key` for keys present in the *current* config, so a directory
whose key was renamed/removed is silently orphaned rather than deleted.

**Impact assessment: harmless but should be cleaned up.** `Snakefile`'s `FAMILIES =
list(data.get("tier1", {}).keys())` (confirmed by direct read in prior audits, re-confirmed here
via the family-key enumeration in §4) does not include the bare `chitinases_gh19` key, so no rule
will read this directory. It contains no data that isn't a strict subset of the current
`chitinases_gh19_class_iv` + `chitinases_gh19_class_i` directories (same 8 accessions, same
sequences — verified by comparing headers). It is not a contamination risk (none of the 4 flagged
bad accessions appear in it either). It is, however, stale disk clutter that could confuse a
future contributor grepping `results/sequences/` for `chitinases_gh19`. Recommend deleting it in
a future pass; not blocking this verification.

---

## Verdict: **CONFIRMED**

All primary claims in the T5 changelog entry hold up under independent re-verification:

1. Zero hits anywhere in `results/sequences/` for all 4 previously-flagged contaminant/wrong-enzyme
   accessions.
2. All 4 spot-checked replacement accessions (`GAV73214.1`, `GAV73783.1`, `XP_019076990.1`,
   `XP_037427783.1`) were independently re-fetched from live NCBI and match the repo's FASTA
   content exactly, residue-for-residue, with organism/annotation consistent with their intended
   role. `Triticum_dicoccoides.fa` is correctly labeled; no stray `Triticum_aestivum.fa` exists.
3. Header format `>{Species}|{acc} family=... species=... acc=... len=...` holds for all 61 FASTA
   files across all 8 families with zero exceptions — not just the family originally flagged.
4. Sequence counts match config non-TODO accession counts for 5/7 spot-checked families exactly;
   the 2 apparent mismatches (`aspartic_proteases_a1b_homology`, `purple_acid_phosphatase`) are
   both fully explained as correct, intentional median-length-filter drops of genuine short
   outlier sequences (confirmed by reproducing the filter arithmetic and, for one, independently
   re-fetching the dropped accession's exact length from NCBI).
5. The `methods_benchmark`/`tier1` merge fix in `fetch_sequences.py` is mechanically correct, has
   no key collisions in the current config, and both `nepenthesins/` and `neprosins/` contain
   fresh, valid, non-empty output.
6. `results/` is fully gitignored and untracked; purging/regenerating it required no special git
   handling.

One non-blocking cleanliness gap noted: a stale, orphaned `results/sequences/chitinases_gh19/`
directory (pre-dating the T4 class split) was not purged by the re-fetch. It is inert (unread by
any Snakemake rule) and contains no contamination, but should be deleted in a future pass.

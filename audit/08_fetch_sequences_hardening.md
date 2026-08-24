# 08 — Harden `fetch_sequences.py` error handling

Code review of `workflow/scripts/fetch_sequences.py` (last touched in commit `96ab469`, which
added `methods_benchmark:` support) found five places where the script would silently discard
data or silently no-op instead of failing loudly. All five are fixed below. Every fix was verified
by actually running the CLI (or the real fetch path) against either the real
`config/enzyme_families.yaml` or a disposable scratch YAML/output dir — not just read-and-reason.

---

## Fix 1 (HIGH): malformed `accessions:` entries vanished with zero log output

**Before:** `if not isinstance(acc_list, list): continue` — a non-list `accessions` entry for a
species (e.g. a YAML edit that turns `Species: [acc1, acc2]` into `Species: acc1` or a nested
mapping) silently dropped that species from the family. No log line, no counter increment, exit 0.

**After** (`workflow/scripts/fetch_sequences.py` around line 216):

```python
for species, acc_list in species_map.items():
    if not isinstance(acc_list, list):
        logger.error(
            "MALFORMED  %s / %s: accessions must be a list, got %s — skipping",
            family_key, species, type(acc_list).__name__,
        )
        n_failed += 1
        continue
```

**Verification command** (throwaway YAML at repo root, `scratch_test_malformed.yaml`, deleted
after the test — never touched `config/enzyme_families.yaml`):

```yaml
tier1:
  test_family:
    expected_length_aa: [100, 200]
    accessions:
      Good_species:
        - BAD07474.1
      Malformed_scalar_species: BAD07474.1
      Malformed_mapping_species:
        nested: BAD07474.1
```

```
$ python3 workflow/scripts/fetch_sequences.py -f scratch_test_malformed.yaml \
    -o /tmp/scratch_fetch_test_malformed -e "$NCBI_EMAIL"
```

**Real output:**

```
2026-08-24 14:18:47,680 INFO __main__: ━━━ Family: test_family (expected 100–200 aa) ━━━
2026-08-24 14:18:49,034 INFO __main__:   FETCHED  Good_species | BAD07474.1 (437 aa)
2026-08-24 14:18:49,034 ERROR __main__: MALFORMED  test_family / Malformed_scalar_species: accessions must be a list, got str — skipping
2026-08-24 14:18:49,034 ERROR __main__: MALFORMED  test_family / Malformed_mapping_species: accessions must be a list, got dict — skipping
2026-08-24 14:18:49,035 INFO __main__: Family test_family: 1 sequences, median 437 aa, length cutoff 350 aa
2026-08-24 14:18:49,055 INFO __main__:   WROTE    Good_species → ...\test_family\Good_species.fa (1 record(s))
2026-08-24 14:18:49,055 INFO __main__: Done. fetched=1  skipped/filtered=0  failed=2
2026-08-24 14:18:49,055 ERROR __main__: 2 accession(s) failed — re-run with --verbose to debug; check NCBI/UniProt availability for those IDs
EXIT CODE: 1
```

Both malformed entries now produce ERROR log lines, increment `n_failed`, and the script exits 1.

---

## Fix 2 (HIGH): `--families` silently accepted nonexistent family keys

**Before:** `family_filter` was never validated against `all_families.keys()`. Passing a stale key
(e.g. `chitinases_gh19`, the pre-`96ab469` name for what is now `chitinases_gh19_class_iv`)
silently produced `fetched=0 skipped/filtered=0 failed=0`, exit 0 — indistinguishable from a
legitimately-empty run.

**After** (around line 181):

```python
if family_filter:
    unknown = family_filter - set(all_families)
    if unknown:
        raise click.BadParameter(
            f"Unknown family key(s): {sorted(unknown)}. Available: {sorted(all_families)}"
        )
    logger.info("Restricting to families: %s", sorted(family_filter))
```

**Verification command** (real config, scratch output dir, deleted after):

```
$ python3 workflow/scripts/fetch_sequences.py -f config/enzyme_families.yaml \
    -o /tmp/scratch_fetch_test -e "$NCBI_EMAIL" --families chitinases_gh19
```

**Real output:**

```
Usage: fetch_sequences.py [OPTIONS]
Try 'fetch_sequences.py --help' for help.

Error: Invalid value: Unknown family key(s): ['chitinases_gh19']. Available: ['aspartic_proteases_a1b_homology', 'chitinases_gh19_class_i', 'chitinases_gh19_class_iv', 'nepenthesins', 'neprosins', 'purple_acid_phosphatase', 'rnase_t2']
EXIT CODE: 2
```

The old key now errors immediately and clearly instead of silently exiting 0.

---

## Fix 3 (HIGH): fetched records were relabeled with the requested accession, destroying evidence of a mismatch

**Before:** `_fetch()`'s caller unconditionally overwrote `record.id`, `record.name`, and
`record.description` with the requested `species`/`acc` string. If NCBI/UniProt resolved a
secondary/redirected accession and returned a *different* protein than requested, all evidence of
that mismatch was silently destroyed — exactly the failure mode that would let a "wrong enzyme
entirely" contamination go undetected (see `audit/05_live_refetch_verification.md`, which had to
catch this class of bug independently after the fact).

**After** (around line 233), a mismatch check runs *before* the record is relabeled:

```python
# Plain GenBank/DDBJ IDs have no pipes (e.g. "BAD07474.1"). Swiss-Prot-sourced
# records — whether from UniProt REST or NCBI's efetch of a Swiss-Prot entry —
# use "sp|ACCESSION|ENTRY_NAME"; the accession is field index 1, NOT the last
# field (verified live: last field is the entry name, e.g. "ASPR_HORVU").
id_parts = record.id.split("|")
fetched_acc = id_parts[1] if len(id_parts) >= 2 else record.id
if fetched_acc.split(".")[0] != acc.split(".")[0]:
    logger.error(
        "ACCESSION MISMATCH %s / %s: requested %s, got %s (%s) — excluding",
        family_key, species, acc, record.id, record.description[:80],
    )
    n_failed += 1
    continue
```

**Important correction found during verification:** the task brief's suggested extraction
(`record.id.split("|")[-1]`) is wrong and would have caused false positives on every legitimate
UniProt fetch. Live-fetching `P42210` (a real UniProt accession in
`config/enzyme_families.yaml`) via the actual UniProt REST endpoint shows:

```
$ python3 -c "... requests.get('https://rest.uniprot.org/uniprotkb/P42210.fasta') ..."
id=           'sp|P42210|ASPR_HORVU'
description=  'sp|P42210|ASPR_HORVU Phytepsin OS=Hordeum vulgare OX=4513 PE=1 SV=1'
```

`split("|")[-1]` gives `'ASPR_HORVU'` (the entry name), not `'P42210'` (the accession) — that
would have incorrectly flagged every Swiss-Prot-format fetch as a mismatch. The same
`sp|ACCESSION|NAME` shape also appears when NCBI's own `efetch` resolves a Swiss-Prot-sourced
accession (confirmed live for `O24603.1`: `id='sp|O24603.1|CHI_ARATH'`), while plain GenBank/DDBJ
accessions have no pipes at all (confirmed live for `BAD07474.1`: `id='BAD07474.1'`). The fix
extracts field index 1 when the id is pipe-delimited, falling back to the whole id otherwise.

**Verification command** (real fetch, real config, scratch output dir):

```
$ python3 workflow/scripts/fetch_sequences.py -f config/enzyme_families.yaml \
    -o /tmp/scratch_fetch_test3 -e "$NCBI_EMAIL" -k "$NCBI_API_KEY" \
    --families nepenthesins --verbose
```

**Real output (tail):**

```
2026-08-24 14:18:24,206 DEBUG __main__: → UniProt: P42210
2026-08-24 14:18:28,752 INFO __main__:   FETCHED  Hordeum_vulgare | P42210 (508 aa)
...
2026-08-24 14:18:30,678 INFO __main__: Family nepenthesins: 14 sequences, median 442 aa, length cutoff 354 aa
...
2026-08-24 14:18:30,691 INFO __main__: Done. fetched=14  skipped/filtered=0  failed=0
EXIT CODE: 0
```

All 14 nepenthesins-family accessions (mix of NCBI GenBank/RefSeq and one UniProt Swiss-Prot
accession, `P42210`) fetched with `failed=0` — confirming the corrected mismatch check does not
false-positive on legitimate fetches, including the exact UniProt-format case that the naive
`[-1]` extraction would have broken.

Credentials were read from `.env` (`NCBI_EMAIL`, `NCBI_API_KEY`) via a small Python snippet that
strips the surrounding quotes/whitespace present in that file; the key value itself was never
printed to any log or terminal output.

---

## Fix 4 (MEDIUM): silent key-collision resolution between `methods_benchmark` and `tier1`

**Before:** `all_families = {**config.get("methods_benchmark", {}), **config.get("tier1", {})}`
with only a code comment ("tier1 wins") — no runtime check. A future accidental duplicate key
across the two sections would silently vanish the `methods_benchmark` definition.

**After** (around line 173):

```python
mb_families = config.get("methods_benchmark", {})
t1_families = config.get("tier1", {})
collision = mb_families.keys() & t1_families.keys()
if collision:
    raise ValueError(
        f"family key(s) defined in both tier1 and methods_benchmark: {sorted(collision)}"
    )
all_families: dict = {**mb_families, **t1_families}
```

**Verification against the real config:**

```
$ python3 -c "
import yaml
config = yaml.safe_load(open('config/enzyme_families.yaml', encoding='utf-8'))
mb, t1 = config.get('methods_benchmark', {}), config.get('tier1', {})
print('methods_benchmark keys:', list(mb.keys()))
print('tier1 keys:', list(t1.keys()))
print('collision:', mb.keys() & t1.keys())
"
methods_benchmark keys: ['nepenthesins', 'neprosins']
tier1 keys: ['chitinases_gh19_class_iv', 'chitinases_gh19_class_i', 'purple_acid_phosphatase', 'rnase_t2', 'aspartic_proteases_a1b_homology']
collision: set()
```

Also confirmed end-to-end via the real CLI against the real config (`--families rnase_t2`
completed normally, exit 0, no `ValueError` raised) — the two sections are genuinely disjoint
today, so this check does not fire; it is a guard against future drift.

---

## Fix 5 (MEDIUM): missing `expected_length_aa` silently disabled all length QC

**Before:** `exp_min, exp_max = family_data.get("expected_length_aa", [0, 99999])` — a family
without this key got `exp_min = 0`, silently making every length-based warning a no-op with no
indication QC was skipped.

**After** (around line 201):

```python
if "expected_length_aa" not in family_data:
    logger.warning(
        "Family %s has no expected_length_aa — length QC disabled for this family",
        family_key,
    )
exp_min, exp_max = family_data.get("expected_length_aa", [0, 99999])
```

**Verification against the real config** — checked every currently-fetched family (all of
`tier1` + `methods_benchmark`):

```
nepenthesins                       HAS [420, 450]
neprosins                          HAS [370, 470]
chitinases_gh19_class_iv           HAS [250, 350]
chitinases_gh19_class_i            HAS [300, 330]
purple_acid_phosphatase            HAS [420, 680]
rnase_t2                           HAS [200, 270]
aspartic_proteases_a1b_homology    HAS [280, 650]
```

**Finding:** none of the 7 currently-fetched families (tier1 + methods_benchmark) lack
`expected_length_aa` — this fix does not currently fire for any real family. It is prophylactic,
guarding against a future `tier2` family (cysteine proteases, TLPs, GH17 glucanases, lipid
transfer proteins, esterases/lipases — none of which are fetched yet per CLAUDE.md Section 3)
being added to `enzyme_families.yaml` without this key.

---

## Summary of execution evidence

| Fix | Real command run | Real result |
|-----|------------------|-------------|
| 1 | CLI against throwaway malformed-accessions YAML, scratch output dir | 2 ERROR log lines, `n_failed=2`, exit 1 |
| 2 | CLI against real config with stale `--families chitinases_gh19`, scratch output dir | `click.BadParameter`, exit 2, clear message |
| 3 | CLI against real config `--families nepenthesins` (real NCBI + UniProt fetch), scratch output dir | 14/14 fetched, `failed=0`; corrected a self-introduced bug in the mismatch-extraction logic before it shipped |
| 4 | Direct config introspection + CLI against real config `--families rnase_t2` | No collision in real config; check does not false-positive |
| 5 | Direct config introspection across all 7 fetched families | None currently missing the key; fix is prophylactic for future `tier2` families |

No scratch files were left behind: `/tmp/scratch_fetch_test*` directories and
`scratch_test_malformed.yaml` were deleted after verification. `results/sequences/` (the real
pipeline output) was never touched by any of these tests.

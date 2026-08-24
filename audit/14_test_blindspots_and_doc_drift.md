# T14: Close test blind spots in test_config_consistency.py; fix CLAUDE.md/README.md tier-list drift

Two independent fixes. Environment: `python3` (3.13) on PATH with `pyyaml`, `click`, `biopython`
installed; `pytest` is NOT installed in this environment (consistent with prior audits, e.g.
audit/10). All verification below was therefore done by manually invoking each test function via
`python3` (importing `tests/test_config_consistency.py` as a module and calling each `test_*`
function directly, matching the manual-verification convention already established in
audit/10_config_gaps_and_consistency_test.md), plus standalone `python3` scripts reproducing each
test's logic against scratch copies of the config files. No real config file was edited for any
scratch test — all scratch files were written under a throwaway `audit/scratch_t14/` directory,
deleted after verification completed.

---

## Fix 1: `tests/test_config_consistency.py` blind spots

### (a) Species listed in both `carnivorous_species:` and `outgroup_species:`

**The gap.** `_known_species_keys()` returns `set(carnivorous_species.keys()) | set(outgroup_species.keys())`
— a plain set union. If a species name appeared as a key in *both* sections, the union still
contains that name exactly once; nothing about a union can distinguish "defined once, correctly"
from "defined twice, silently overwritten." Meanwhile `detect_convergence.py`'s
`_load_species_meta()` (workflow/scripts/detect_convergence.py:37-85) loads `carnivorous_species`
first, then unconditionally overwrites `meta[species_name]` for each `outgroup_species` entry — so
a species in both sections ends up classified as non-carnivorous with `carnivory_origin=0`, silently
dropped from every convergence count.

**Verification steps performed:**
1. Built a scratch copy of `config/species.yaml` with `Cephalotus_follicularis`'s real entry
   duplicated verbatim into `outgroup_species:` (real file untouched).
2. Ran the logic of all 5 pre-existing tests against this scratch file. Result: **all 5 passed**
   — none detected the duplication. Directly inspected
   `_load_species_meta(scratch_path)["Cephalotus_follicularis"]` and confirmed it resolved to
   `{'carnivorous': False, 'carnivory_origin': 0, 'name': 'Cephalotus_follicularis'}` — i.e. the
   bug is real and silent, exactly as hypothesized.
3. Added `test_carnivorous_and_outgroup_species_disjoint()`, asserting
   `set(carnivorous_species.keys()) & set(outgroup_species.keys()) == set()`.
4. Ran the new assertion's logic against the scratch duplicate file: `overlap = {'Cephalotus_follicularis'}`
   → **fails**, correctly catching the injected bug.
5. Ran it against the real `config/species.yaml`: `overlap = set()` → **passes**. No duplication
   exists today.

### (b) Duplicate `code:` values across both sections

**The gap.** `_load_species_meta()` also registers each species under `meta[code]`
(detect_convergence.py:63-65, 75-77), in the same carnivorous-then-outgroup load order. Two
species sharing a `code` would silently overwrite the same way as (a).

**Verification steps performed:**
1. Confirmed via a direct `Counter` over every `code:` value in the real `config/species.yaml`
   (both sections combined) that there are currently **zero duplicates**.
2. Built a scratch copy of `species.yaml` with `Drosera_binata`'s `code` changed from `Dbin` to
   `Dade` (colliding with `Drosera_adelae`'s real code).
3. Added `test_species_codes_are_unique()`, which walks both sections, tracks `code -> first
   species name seen`, and collects any second-or-later species sharing that code into a `dupes`
   dict.
4. Ran the new assertion's logic against the scratch dup-code file: `dupes = {'Dade': ['Drosera_adelae', 'Drosera_binata']}`
   → **fails**, correctly catching the injected bug.
5. Ran it against the real file: `dupes = {}` → **passes**.

### (c) Every `tier1:` family spans `>= convergence.min_lineages` distinct carnivory origins (non-TODO accessions)

**The gap.** No existing test enforced the exact invariant that motivated moving
`nepenthesins`/`neprosins` (2026-08-22) and `chitinases_gh19_class_i` (2026-08-24) out of
`tier1:` — a `tier1:` family whose real (non-TODO) accessions come from fewer than
`config.yaml`'s `convergence.min_lineages` distinct `carnivory_origin` values can mathematically
never produce a cross-lineage convergence result (`detect_convergence.py`'s
`_detect_convergent_sites()` explicitly does `if n_origins < min_lineages: continue`). This can
silently regress if someone edits accessions later (e.g. TODO-ing out the last real accession
from a second origin) without noticing the family no longer qualifies for `tier1:`.

**Verification steps performed:**
1. Read the real `convergence.min_lineages` value from `config/config.yaml`: **2** (did not
   hardcode this — the new test reads it live via `config["convergence"]["min_lineages"]`).
2. Wrote a standalone script computing, for each real `tier1:` family, the set of distinct
   `carnivory_origin` values among species whose accessions list contains at least one non-`TODO`
   entry (species not in `carnivorous_species:` — i.e. outgroups — are excluded from the origin
   count, since only carnivorous-lineage origins count toward cross-lineage convergence). Result
   against the real config:
   - `chitinases_gh19_class_iv`: origins `{1, 2, 3}` (n=3) → PASS
   - `purple_acid_phosphatase`: origins `{1, 2, 3}` (n=3) → PASS
   - `rnase_t2`: origins `{1, 2}` (n=2) → PASS (exactly at the threshold)
   - `aspartic_proteases_a1b_homology`: origins `{1, 2, 3}` (n=3) → PASS
   All 4 real `tier1:` families currently pass.
3. Added `test_tier1_families_span_min_lineages_carnivory_origins()` implementing this exact
   logic (with a `_is_todo_only()` helper matching the accession-list TODO-sentinel convention
   used throughout `config/enzyme_families.yaml`).
4. Built a scratch copy of `enzyme_families.yaml` with `aspartic_proteases_a1b_homology`'s
   `Cephalotus_follicularis` and `Sarracenia_purpurea` accession entries replaced with `[TODO]`
   (simulating someone TODO-ing out the last real accessions for origins 2 and 3, leaving only
   origin 1 — Drosera/Dionaea).
5. Ran the new test's logic against the scratch file: `aspartic_proteases_a1b_homology` origins
   `{1}` (n=1) < min_lineages=2 → **fails**, correctly catching the injected regression. The other
   3 families (untouched in the scratch copy) still pass.
6. Ran it against the real, unmodified `config/enzyme_families.yaml`: all 4 families pass.

### (d) `test_seq_id_to_species_roundtrip_on_species_names` docstring overclaim, fixed

**The gap.** The test's docstring claims it guards the "species-name-vs-code keying bug" in
`_seq_id_to_species()`/`_load_species_meta()`, but its body built `known_species_names` via
`_known_species_keys(species)` — reading `config/species.yaml` directly and independently of
`module._load_species_meta()`'s actual output. It therefore only ever exercised
`_seq_id_to_species()` in isolation with a hand-built, always-correct `known_species` set; it
never actually called `_load_species_meta()` and could not detect a regression in the loader
itself. (`test_load_species_meta_keys_include_full_species_names`, a sibling test, does correctly
guard the loader directly — so this was a docstring-overclaim / test-scope mismatch, not an
uncovered loader regression, but it's still misleading and was fixed per the task.)

**Fix applied:** changed the test to derive its `known_species` argument from
`set(module._load_species_meta(SPECIES_YAML).keys())` (the loader's actual output), while still
iterating the round-trip subject names from `_known_species_keys(species)` (the real species
names read directly from `species.yaml`, used only to build each fake `"{name}|FAKEACC.1"` seq
ID — this list of names-to-test is independent of the bug under test).

**Verification steps performed:**
1. Built a scratch copy of `workflow/scripts/detect_convergence.py` with `_load_species_meta()`'s
   name-key registration lines removed — i.e. reverted to registering each species **only** under
   its short `code`, not its full species name (simulating "before the species-keying fix" state
   described in `audit/07_detect_convergence_deep_fix.md`).
2. Ran the **pre-fix** (docstring-overclaiming) test logic against this reverted module: **0
   failures reported → test PASSES**, even though the loader is broken. This directly confirms
   the overclaim — the old test's body never actually calls the reverted `_load_species_meta()`,
   so it cannot see the bug.
3. Ran the **fixed** test logic (using `module._load_species_meta(SPECIES_YAML).keys()`) against
   the same reverted module: **29 failures** (one per species, e.g.
   `('Aldrovanda_vesiculosa|FAKEACC.1', None)` — every lookup fails because the reverted loader's
   key set contains only 4-letter codes, none of which prefix-match a
   `"Full_Species_Name|FAKEACC.1"` string) → **test correctly FAILS**, exactly as required.
4. Ran the fixed test logic against the real, unmodified `detect_convergence.py`: **0 failures →
   PASSES**.

### Full test-suite run against the real config

All 8 test functions (5 original — one fixed, one unchanged aside from the docstring fix — plus 3
new) were executed directly via `python3` (since `pytest` isn't installed here) by importing
`tests/test_config_consistency.py` as a module and calling each `test_*` function in turn:

```
Found 8 test functions:
  PASS  test_all_accession_species_exist_in_species_yaml
  PASS  test_carnivorous_and_outgroup_species_disjoint
  PASS  test_carnivorous_species_have_nonzero_carnivory_origin
  PASS  test_load_species_meta_keys_include_full_species_names
  PASS  test_seq_id_to_species_roundtrip_on_species_names
  PASS  test_species_codes_are_unique
  PASS  test_tier1_and_methods_benchmark_families_disjoint
  PASS  test_tier1_families_span_min_lineages_carnivory_origins
```

`python3 -m py_compile tests/test_config_consistency.py` confirmed the file is syntactically valid.

---

## Fix 2: CLAUDE.md / README.md tier-list drift

**The gap.** `config/enzyme_families.yaml`'s actual keys (verified via live PyYAML parse):

```
tier1:             ['chitinases_gh19_class_iv', 'purple_acid_phosphatase', 'rnase_t2', 'aspartic_proteases_a1b_homology']
methods_benchmark: ['chitinases_gh19_class_i', 'nepenthesins', 'neprosins']
```

`CLAUDE.md` §3 still listed `chitinases_gh19_class_i` as Tier 1 item #2 (5-item numbered list),
directly violating its own stated rule ("This list must match `config/enzyme_families.yaml`'s
actual `tier1:` keys exactly — do not let it drift") left in place from the 2026-08-22 restructure
but never updated when `chitinases_gh19_class_i` moved to `methods_benchmark:` on 2026-08-24
(audit/10_config_gaps_and_consistency_test.md). `README.md`'s family table had the same gap: its
GH19 Class I row used the bare key `chitinases_gh19_class_i` while its sibling
`methods_benchmark` rows (nepenthesins, neprosins) correctly show the `methods_benchmark.`
prefix; the lead-in paragraph above the table also still said "GH19 Class I is retained tier1."

**Fixes applied:**
- `CLAUDE.md` §3: renumbered Tier 1 to the real 4 items (`chitinases_gh19_class_iv`,
  `purple_acid_phosphatase`, `rnase_t2`, `aspartic_proteases_a1b_homology`); moved
  `chitinases_gh19_class_i` into the "Methods Benchmark" subsection as its own bullet, using the
  same single-origin/`min_lineages` framing already used for `nepenthesins`/`neprosins`, plus its
  own domain-architecture note (retained from the original entry); generalized the Methods
  Benchmark intro blurb from "Nepenthes-only families" to "Nepenthes-only, and
  single-Caryophyllales-origin, families" since GH19 Class I is single-origin but not
  Nepenthes-only (it's 3 *Drosera* species).
- `README.md`: changed the GH19 Class I table row's key to
  `` `methods_benchmark.chitinases_gh19_class_i` `` (matching the nepenthesins/neprosins row
  convention) and its "Convergence target?" cell to "No — single-origin structural comparison
  only" (matching the sibling rows' "No — single-origin methods benchmark" framing); rewrote the
  lead-in blockquote above the table, which previously said "GH19 Class I is retained tier1,"
  to state all three (GH19 Class I, nepenthesins, neprosins) live under `methods_benchmark:` and
  explain each one's single-origin reason, and added a note that GH19 Class I's move happened
  2026-08-24 per audit/10.

**Verification (programmatic, not eyeballed):** wrote a `python3` script that (1) PyYAML-parses
`config/enzyme_families.yaml` for the real `tier1`/`methods_benchmark` key sets; (2) extracts
CLAUDE.md's "Tier 1" and "Methods Benchmark" subsections by regex and checks which of the 7 total
family keys appear in backtick-quoted form (`` `key` ``) in each; (3) checks README.md for
backtick-quoted bare keys and `methods_benchmark.`-prefixed keys. Results:

```
config tier1: ['aspartic_proteases_a1b_homology', 'chitinases_gh19_class_iv', 'purple_acid_phosphatase', 'rnase_t2']
config methods_benchmark: ['chitinases_gh19_class_i', 'nepenthesins', 'neprosins']

CLAUDE.md Tier 1 section backtick-keys found: ['aspartic_proteases_a1b_homology', 'chitinases_gh19_class_iv', 'purple_acid_phosphatase', 'rnase_t2']
  matches config tier1? True
CLAUDE.md Methods Benchmark section backtick-keys found: ['chitinases_gh19_class_i', 'nepenthesins', 'neprosins']
  matches config methods_benchmark? True

README.md has 'methods_benchmark.chitinases_gh19_class_i' prefix: True
README.md has 'methods_benchmark.nepenthesins' prefix: True
README.md has 'methods_benchmark.neprosins' prefix: True
README.md tier1 keys incorrectly prefixed methods_benchmark.: False (all 4 checked)
```

Both documents now exactly match `config/enzyme_families.yaml`'s real `tier1:`/`methods_benchmark:`
partition.

---

## Files changed

- `tests/test_config_consistency.py` — fixed test (d)'s docstring/body mismatch; added
  `_load_config()` helper and 3 new test functions ((a), (b), (c) above); updated module
  docstring to list all guarded bugs (a)-(g).
- `CLAUDE.md` — §3 Tier 1 list renumbered to 4 items; `chitinases_gh19_class_i` moved to Methods
  Benchmark subsection.
- `README.md` — GH19 Class I table row and lead-in paragraph corrected to reflect
  `methods_benchmark:` membership.

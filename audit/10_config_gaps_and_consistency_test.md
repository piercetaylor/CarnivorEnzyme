# 10 — Fix config structural gaps, add consistency test

Task: (1) `chitinases_gh19_class_i` was left in `tier1:` despite being explicitly
documented as structural-comparison-only, guaranteeing a header-only empty
convergence-detection output on every run; (2) 6 species used as `accessions:` keys in
`config/enzyme_families.yaml` have no matching `config/species.yaml` entry, plus one
malformed species key; (3) write the `tests/test_config_consistency.py` regression test a
prior review recommended. All fixes applied to the real config files; every assertion
verified by actually running it (`python3`, real PyYAML parsing) against the real files —
not just read-and-reasoned.

---

## Fix 1 — `chitinases_gh19_class_i` moved out of `tier1:`

**Verified the premise directly** before touching anything:

```bash
python3 -c "
import yaml
fam = yaml.safe_load(open('config/enzyme_families.yaml', encoding='utf-8'))
sp = yaml.safe_load(open('config/species.yaml', encoding='utf-8'))
species = fam['tier1']['chitinases_gh19_class_i']['accessions'].keys()
for s in species:
    print(s, '->', sp['carnivorous_species'][s]['carnivory_origin'])
"
```

Output (pre-fix): `Drosera_capensis -> 1`, `Drosera_adelae -> 1`, `Drosera_binata -> 1`
— all three species are `carnivory_origin: 1` (Caryophyllales, per the single-origin
scheme; `config/species.yaml` header, `audit/02_carnivory_origin_reassignment.md`).
`config/config.yaml`'s `convergence.min_lineages: 2` means a family entirely within one
origin mathematically cannot pass the `n_origins >= min_lineages` check in
`detect_convergence.py`'s `_detect_convergent_sites()` — every run against this family
produces a header-only TSV, exit 0, no error. `Snakefile:44-49` derives
`FAMILIES = list(data.get("tier1", {}).keys())`, which feeds `phase2`, `phase3a`,
`phase4a`, `phase5d`, and `rule all` — so this family was silently baked into every one of
those targets.

**Precedent**: `nepenthesins` and `neprosins` were already moved from `tier1:` to a new
`methods_benchmark:` section on 2026-08-22 for the exact same single-origin reason (see
`audit/03_merops_restructure_and_neprosin_rescope.md`). This fix applies the identical
treatment to `chitinases_gh19_class_i` for consistency — it was already labeled
`note: "Structural comparison only — not used for cross-lineage convergence detection"` in
its own YAML but the label had never been backed up by actually moving the entry.

**Applied**: moved the entire `chitinases_gh19_class_i` block (all comments, both
top-of-section rationale comments and per-accession comments, all 3 species' accessions)
verbatim from `tier1:` into `methods_benchmark:`, placed at the top of that section (before
`nepenthesins`). Left a placeholder comment in `tier1:` at the old location cross-referencing
the move, mirroring the style of the existing "4. NEPENTHESINS — MOVED" placeholder left
behind by the 2026-08-22 move. No accessions, TODOs, or notes were dropped.

**Verification**:

```bash
python3 -c "
import yaml
d = yaml.safe_load(open('config/enzyme_families.yaml', encoding='utf-8'))
print('tier1:', sorted(d['tier1'].keys()))
print('methods_benchmark:', sorted(d['methods_benchmark'].keys()))
"
```

Output (post-fix):
```
tier1: ['aspartic_proteases_a1b_homology', 'chitinases_gh19_class_iv', 'purple_acid_phosphatase', 'rnase_t2']
methods_benchmark: ['chitinases_gh19_class_i', 'nepenthesins', 'neprosins']
```

`chitinases_gh19_class_i` is now under `methods_benchmark`, absent from `tier1` — matches
the requested check exactly.

**Downstream reference check**: grepped the whole repo for `chitinases_gh19_class_i\b`.
Found 8 files; only `config/enzyme_families.yaml` (fixed above) and
`workflow/scripts/braker2/04_run_hmmer_scan.sh` are non-documentation/non-audit hits. The
`04_run_hmmer_scan.sh` hit is a comment only (`... being filed under
chitinases_gh19_class_iv vs chitinases_gh19_class_i — this script only finds candidates,
it doesn't classify them`), not a `FAMILY_HMMS`/`MIN_LEN` array key — that script's actual
family-keyed arrays only reference `chitinases_gh19_class_iv`, `purple_acid_phosphatase`,
`rnase_t2`, `aspartic_proteases_a1b_homology`, so no fix needed there. The remaining 6 hits
are narrative audit/README/CLAUDE.md documentation of the historical restructuring, out of
scope.

**`config/config.yaml` / `config/substrates.yaml` check**: `chitinases_gh19_class_i` was
never a docking or FEP target (it's structural-comparison-only — no catalytic-site
convergence data to dock against or run FEP on), so its correct state is *absence* from
`fep.substrates_for_fep` (config.yaml) and `substrates.yaml`. Confirmed both files only
reference `chitinases_gh19_class_iv` (already correctly updated in the 2026-08-22 GH19
class-split pass, `audit/04_gh19_class_split.md`) and have no `chitinases_gh19_class_i`
entry — this is the correct, intentional state, not a gap.

---

## Fix 2 — species.yaml gaps + malformed key

**Gap-finding script** (written to `audit/_tmp_species_gap_check.py`, run, then deleted —
its logic is now permanently encoded as
`tests/test_config_consistency.py::test_all_accession_species_exist_in_species_yaml`):

```python
import yaml
from pathlib import Path

ROOT = Path(__file__).parent.parent
fam = yaml.safe_load(open(ROOT / "config/enzyme_families.yaml", encoding="utf-8"))
spec = yaml.safe_load(open(ROOT / "config/species.yaml", encoding="utf-8"))

used_species = set()
for section in ("tier1", "methods_benchmark"):
    for family_name, family_data in fam.get(section, {}).items():
        accessions = family_data.get("accessions", {}) or {}
        for sp in accessions.keys():
            used_species.add(sp)

known_species = set(spec.get("carnivorous_species", {}).keys()) | set(spec.get("outgroup_species", {}).keys())

missing = sorted(used_species - known_species)
print("Species referenced in enzyme_families.yaml accessions but missing from species.yaml:")
for sp in missing:
    print(" -", sp)
print()
print(f"Total used species keys: {len(used_species)}")
print(f"Total known species.yaml keys: {len(known_species)}")
print(f"Missing count: {len(missing)}")
```

**Output — BEFORE fix**:
```
Species referenced in enzyme_families.yaml accessions but missing from species.yaml:
 - Drosera_rotundifolia
 - Nepenthes_ampullaria
 - Nepenthes_bicalcarata
 - Nepenthes_mirabilis_PAP
 - Nepenthes_rafflesiana
 - Nepenthes_ventrata
 - Triticum_dicoccoides

Total used species keys: 29
Total known species.yaml keys: 23
Missing count: 7
```

Matches the task's predicted list exactly (7 gaps: 6 real missing species + 1 malformed
key).

**Output — AFTER fix**:
```
Species referenced in enzyme_families.yaml accessions but missing from species.yaml:

Total used species keys: 28
Total known species.yaml keys: 29
Missing count: 0
```

(Used-key count dropped 29→28 because `Nepenthes_mirabilis_PAP` was renamed to
`Nepenthes_mirabilis`, which was already a used key elsewhere, collapsing the union by
one; known-key count rose 23→29 from the 6 new `species.yaml` entries.)

### Per-species reasoning

**`Nepenthes_mirabilis_PAP` (malformed key, not a real gap)** — `tier1.purple_acid_phosphatase`
had a species-key entry literally named `Nepenthes_mirabilis_PAP` (a species name with a
family-disambiguation suffix awkwardly appended — not a valid species name and not how any
other family in this file keys its accessions). Checked first whether a plain
`Nepenthes_mirabilis` key already existed elsewhere in `purple_acid_phosphatase`'s own
`accessions:` mapping (it did not — grepped all 4 occurrences of
`Nepenthes_mirabilis`/`Nepenthes_mirabilis_PAP` in the file; the other 3 are in
*different* families' own `accessions:` dicts — `chitinases_gh19_class_iv`,
`methods_benchmark.nepenthesins`, `methods_benchmark.neprosins` — which is fine, each
family has its own independent accessions mapping). Renamed the key to `Nepenthes_mirabilis`.
This was a plain rename, not a merge (no duplicate-key YAML clobbering risk), since no
sibling key existed in that same dict.

**`Drosera_rotundifolia`** — has real (non-TODO) accessions in `tier1.chitinases_gh19_class_iv`
(`GAB2242148.1`, `GAB2242147.1`) and `tier1.aspartic_proteases_a1b_homology`
(`GAB2239788.1`). Added to `carnivorous_species:` with `carnivory_origin: 1` (same as all
other Caryophyllales species under the single-origin scheme). `taxid: 173423` — looked up
via WebSearch/NCBI Taxonomy Browser
(https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=173423), not invented. For
`genome:`, WebSearch alone did not surface a public GCA/GCF accession for the "Droror1"
gene-ID prefix seen in the accession comments (`Droror1_Dr00018923` etc.) — rather than
guess, fetched the actual NCBI GenBank record for `GAB2242148.1` directly
(`efetch.fcgi?db=protein&id=GAB2242148.1&rettype=gb`) and read its DBSOURCE line, which
gave BioProject `PRJDB18922`, assembly name `Droro_R1_v2`, WGS prefix `BAAGII01`
(139x-coverage Illumina HiSeq X Ten + Nanopore GridION hybrid assembly). Used the
BioProject accession as the `genome:` value (the same convention already used for
`Aldrovanda_vesiculosa: genome: PRJNA918534`, another entry with no GCA/GCF on file) rather
than leaving it blank or guessing a GCA number that was never directly confirmed.
`tier: 1` (has real accessions used by two tier1 families, same treatment as Dionaea/Cephalotus).

**`Triticum_dicoccoides`** — real outgroup accession in `tier1.aspartic_proteases_a1b_homology`
(`XP_037427783.1`) and `methods_benchmark.nepenthesins` (same accession, duplicated per that
section's existing convention of duplicating A1B outgroups across both split families). Added
to `outgroup_species:`, matching the field set used by the existing `Hordeum_vulgare` /
`Secale_cereale` entries (`taxid`, `code`, `genome`, `pitcher_pH: null`, `clade`, `order`,
`trap: null`, `carnivorous: false`, `expression_sra: null`, `tier`, optional `notes`).
`taxid: 85692` and `genome: GCA_002162155.1` (WEWSeq v1.0) — both looked up via WebSearch
against NCBI (taxonomy browser id 85692; NCBI Datasets genome page for
`GCA_002162155.1`), not invented. `genome_source` field omitted to match the existing
outgroup entries' field set exactly (none of them carry a `genome_source:` field — only
`carnivorous_species:` entries do); cited Avni et al. 2017, *Science* 357:93-97 (the WEWSeq
v1.0 genome paper) in the `notes:` field instead.

**`Nepenthes_ampullaria`** — has real, verified accessions (`ARA95696.1`, `ARA95695.1`) in
`methods_benchmark.neprosins`, explicitly called out in that file's own comment as "high-pH
pitcher (pH~5) — useful for pH comparison" for the within-*Nepenthes* case study described
in CLAUDE.md Phase 4A / methods_benchmark rationale. Added with `carnivory_origin: 1`,
`pitcher_pH: 5.0` (matching that in-file comment), `tier: 1`. `taxid: 150949` verified via
NCBI Taxonomy Browser. No genome found via WebSearch as of 2026-08-24 (only RNA-seq data,
PMC4707257) — `genome: none`, consistent with the file's existing convention for species
with transcriptome-only data.

**`Nepenthes_ventrata`** — has a real, verified accession (`C0HLV2.1`) in
`methods_benchmark.neprosins`, explicitly the sequence used for the 7ZVA/7ZVB/7ZVC/7ZU8
crystal structures — the single experimental structure this whole project benchmarks
against (CLAUDE.md §2, §3). This is arguably the single most important species missing
from `species.yaml`. `taxid: 1744888` — *N.* × *ventrata* is a horticultural hybrid name
without its own clean NCBI Taxonomy Browser hit via plain WebSearch (multiple queries
returned only related-species IDs); resolved directly against the NCBI Taxonomy database
itself via `eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=Nepenthes%20x%20ventrata`
(returned a single hit, id `1744888`), then confirmed via `esummary.fcgi?db=taxonomy&id=1744888`
that the scientific name on record is "Nepenthes ventricosa x Nepenthes alata" — the
recognized parentage of the hybrid commonly sold/cited as *N.* × *ventrata*. `carnivory_origin: 1`,
`tier: 1` (crystal-structure-anchor species). Only the mitochondrial genome is published
(PMC7800874) — `genome: none` for the nuclear assembly, noted explicitly.

**`Nepenthes_bicalcarata`, `Nepenthes_rafflesiana`** — both are TODO-only in
`methods_benchmark.neprosins` (no real accession, per that file's own TODO comments: no
public genome/accession for *N. bicalcarata*; *N. rafflesiana* has no NCBI accession and
`methods_benchmark.neprosins` already notes *N.* × *ventrata* (`C0HLV2`) is used as a
surrogate since *N.* × *ventrata* is (in one common parentage account) partly derived from
*N. rafflesiana*). Per the task's own suggested judgment call, gave these two a
**lighter-weight placeholder** rather than a fully-fleshed-out entry: real, WebSearch-verified
`taxid` (330715 and 150990 respectively, both confirmed via NCBI Taxonomy Browser), `genome:
none`, `tier: 2` (vs. `tier: 1` for the species with real accessions), and a `notes:` field
that explicitly flags TODO-only status so a future pass doesn't mistake these for
fully-verified tier1-grade entries.

---

## Fix 3 — `tests/test_config_consistency.py`

Checked `tests/conftest.py` first: it defines `project_root` and `test_data` pytest
fixtures, unused by this test file (no test data files are needed — everything reads the
real `config/*.yaml` and imports the real `workflow/scripts/detect_convergence.py`
directly). Followed the existing `tests/test_*.py` style (plain `def test_...():` functions,
plain `assert` statements) since none of the other test files use fixtures either — they're
all still `TODO: implement` placeholders.

Implemented 5 test functions (the requested (a)–(d), plus one extra direct regression guard
for (d)'s underlying bug):

- `test_all_accession_species_exist_in_species_yaml` — (a)
- `test_carnivorous_species_have_nonzero_carnivory_origin` — (b)
- `test_tier1_and_methods_benchmark_families_disjoint` — (c)
- `test_seq_id_to_species_roundtrip_on_species_names` — (d)
- `test_load_species_meta_keys_include_full_species_names` — extra: directly checks that
  `_load_species_meta()`'s returned dict is keyed by full species names (not just short
  `code`s), which is the actual mechanism (d) depends on.

**Import mechanics note**: `workflow/scripts/detect_convergence.py` is not a package (no
`__init__.py`, not installed) and does `import ete3` at module level. This test
environment's `python3` has `pyyaml`, `click`, `biopython`, `scipy`, `numpy` installed but
**not** `ete3` (confirmed: `python3 -c "import ete3"` → `ModuleNotFoundError`). Since
`_seq_id_to_species()` and `_load_species_meta()` are pure-Python and never touch `ete3`,
the test file loads `detect_convergence.py` by file path via `importlib.util.spec_from_file_location`
and, if the real `ete3` isn't importable, installs a minimal stub `ete3` module into
`sys.modules` first (just enough for the module-level `import ete3` line to succeed) rather
than skipping the test outright. This is test-harness plumbing only — it does not change
or mock anything about the two functions under test.

### Manual verification (pytest not installed in this environment; ran the assertions directly)

```bash
python3 -c "
import sys
sys.path.insert(0, 'tests')
import test_config_consistency as t
t.test_all_accession_species_exist_in_species_yaml()
print('(a) PASS')
t.test_carnivorous_species_have_nonzero_carnivory_origin()
print('(b) PASS')
t.test_tier1_and_methods_benchmark_families_disjoint()
print('(c) PASS')
"
```
Output:
```
(a) PASS
(b) PASS
(c) PASS
```

```bash
python3 -c "
import sys
sys.path.insert(0, 'tests')
import test_config_consistency as t
t.test_seq_id_to_species_roundtrip_on_species_names()
print('(d) PASS')
t.test_load_species_meta_keys_include_full_species_names()
print('(d2) PASS')
"
```
Output:
```
(d) PASS
(d2) PASS
```

All 5 assertions currently pass against the real, fixed config files. (d)/(d2) also
directly confirm that the species-name-vs-code keying fix already landed in
`detect_convergence.py._load_species_meta()` (per `audit/07_detect_convergence_deep_fix.md`,
finding (1)) is still in place and functioning: every one of the 29 species.yaml names
round-trips correctly through `_seq_id_to_species()` when constructed as the
`"{Species_name}|{accession}"` format `fetch_sequences.py` actually writes.

### Full config re-verification after all fixes

```bash
python3 -c "
import yaml
fam = yaml.safe_load(open('config/enzyme_families.yaml', encoding='utf-8'))
sp = yaml.safe_load(open('config/species.yaml', encoding='utf-8'))
print('tier1:', sorted(fam['tier1'].keys()))
print('methods_benchmark:', sorted(fam['methods_benchmark'].keys()))
print('carnivorous_species count:', len(sp['carnivorous_species']))
print('outgroup_species count:', len(sp['outgroup_species']))
"
```
Output:
```
tier1: ['aspartic_proteases_a1b_homology', 'chitinases_gh19_class_iv', 'purple_acid_phosphatase', 'rnase_t2']
methods_benchmark: ['chitinases_gh19_class_i', 'nepenthesins', 'neprosins']
carnivorous_species count: 23
outgroup_species count: 6
```

(23 = 18 pre-existing carnivorous_species entries + 5 new: `Drosera_rotundifolia`,
`Nepenthes_ampullaria`, `Nepenthes_bicalcarata`, `Nepenthes_rafflesiana`,
`Nepenthes_ventrata`. 6 = 5 pre-existing outgroup_species entries + `Triticum_dicoccoides`.)

# Audit: Fix `detect_convergence.py` Species-Metadata Loader

## Summary

Fixed a critical bug in `workflow/scripts/detect_convergence.py` where the `_load_species_meta()` function was reading from a non-existent `species:` key in the YAML, causing it to always load zero species and silently fail to detect convergence.

## The Bug

**File:** `workflow/scripts/detect_convergence.py`  
**Function:** `_load_species_meta()` (lines 37–49, before fix)

**Original code:**
```python
def _load_species_meta(species_yaml: Path) -> dict:
    """Return {species_code: {carnivorous, carnivory_origin, ...}}."""
    with open(species_yaml) as fh:
        data = yaml.safe_load(fh)
    meta = {}
    for entry in data.get("species", []):  # <-- BUG: "species" key does not exist
        code = entry.get("code") or entry.get("name", "").replace(" ", "_")
        meta[code] = {
            "carnivorous": bool(entry.get("carnivorous", False)),
            "carnivory_origin": entry.get("carnivory_origin", 0),
            "name": entry.get("name", code),
        }
    return meta
```

**Root cause:**
- `data.get("species", [])` returns `[]` because `config/species.yaml` has NO top-level `species:` key
- The YAML structure actually has two sections: `carnivorous_species:` and `outgroup_species:`, each containing a dict of species keyed by display name
- The loop iterates zero times, returning an empty dict `{}`
- Downstream code (`_get_carnivorous_branch_lengths()` at line 216, `_seq_id_to_species()` check) silently fails because `species_meta` is empty
- No error is raised — the convergence detection simply finds nothing

**Impact:**
- Convergent substitution detection is completely broken — always returns zero results
- Pipeline produces empty TSV output files with no warning
- The bug was silent and would not be caught without running actual convergence tests

## The Fix

**New code (lines 37–59):**
```python
def _load_species_meta(species_yaml: Path) -> dict:
    """Return {species_code: {carnivorous, carnivory_origin, name}}."""
    with open(species_yaml) as fh:
        data = yaml.safe_load(fh)
    meta = {}

    # Load carnivorous species
    for species_name, entry in data.get("carnivorous_species", {}).items():
        code = entry.get("code") or species_name.replace(" ", "_")
        meta[code] = {
            "carnivorous": True,
            "carnivory_origin": entry.get("carnivory_origin", 0),
            "name": species_name,
        }

    # Load outgroup (non-carnivorous) species
    for species_name, entry in data.get("outgroup_species", {}).items():
        code = entry.get("code") or species_name.replace(" ", "_")
        meta[code] = {
            "carnivorous": False,
            "carnivory_origin": 0,
            "name": species_name,
        }

    return meta
```

**Key changes:**
1. Iterate `data.get("carnivorous_species", {}).items()` to extract species name and metadata dict
2. For each carnivorous species, set `"carnivorous": True` (inferred from section membership, not from an explicit field)
3. Extract `code` from `entry.get("code")`, falling back to derived name if missing
4. Extract `carnivory_origin` from the entry (present only in carnivorous_species)
5. Separately iterate `data.get("outgroup_species", {}).items()`
6. For outgroups, set `"carnivorous": False` and `"carnivory_origin": 0` (not derived, hard-coded)
7. Both sections are merged into the same output dict, maintaining the same return type: `{species_code: {...}}`

**Pattern rationale:**
- Follows the same correct pattern used in `workflow/scripts/root_tree.py` line 31: `config.get("outgroup_species", {}).keys()`
- Infers `carnivorous` from section membership (reliable) rather than checking an explicit field (the YAML does not consistently set this field — carnivorous_species entries do not have `carnivorous: true`, they are identified by being in the section)

## Manual Trace Against Actual YAML

Traced the fix against the complete content of `config/species.yaml` (lines 1–348):

### Carnivorous species section (lines 28–277)
Expected: 18 entries, all with `carnivorous=True`
```
Nepenthes_gracilis     → code=Ngra,  carnivory_origin=1
Nepenthes_alata        → code=Nala,  carnivory_origin=1
Nepenthes_mirabilis    → code=Nmir,  carnivory_origin=1
Drosera_capensis       → code=Dcap,  carnivory_origin=2
Drosera_adelae         → code=Dade,  carnivory_origin=2
Drosera_spatulata      → code=Dspa,  carnivory_origin=2
Drosera_binata         → code=Dbin,  carnivory_origin=2
Dionaea_muscipula      → code=Dmus,  carnivory_origin=3
Aldrovanda_vesiculosa  → code=Aves,  carnivory_origin=3
Cephalotus_follicularis→ code=Cfol,  carnivory_origin=4
Sarracenia_purpurea    → code=Spur,  carnivory_origin=5
Sarracenia_alata       → code=Sala,  carnivory_origin=5
Darlingtonia_californica→code=Dcal,  carnivory_origin=5
Heliamphora_ciliata    → code=Hcil,  carnivory_origin=5
Utricularia_gibba      → code=Ugib,  carnivory_origin=6
Pinguicula_moranensis  → code=Pmor,  carnivory_origin=6
Pinguicula_vulgaris    → code=Pvul,  carnivory_origin=6
Byblis_filifolia       → code=Bfil,  carnivory_origin=7
```
Trace result: **All 18 entries loaded with `carnivorous=True` ✓**

### Outgroup species section (lines 282–348)
Expected: 5 entries, all with `carnivorous=False`
```
Arabidopsis_thaliana → code=Atha, carnivorous=False
Vitis_vinifera       → code=Vvin, carnivorous=False
Solanum_lycopersicum → code=Slyo, carnivorous=False
Hordeum_vulgare      → code=Hvul, carnivorous=False
Secale_cereale       → code=Scer, carnivorous=False
```
Trace result: **All 5 entries loaded with `carnivorous=False` ✓**

### Total dict size
Expected: 23 entries (18 + 5)  
**Trace result: Fixed loader will produce exactly 23 entries ✓**

### Code extraction verification
All entries in both sections have explicit `code:` field (no fallback to derived names needed), so the fallback logic in `code = entry.get("code") or species_name.replace(" ", "_")` is not exercised but is defensive and correct.

## Dead Code Audit

Searched `workflow/scripts/detect_convergence.py` for duplicate or unused species-related loader functions:
- No `_load_family_outgroups()` found
- No other `_load_species*()` functions found
- No dead code to remove in this file

(Note: `root_tree.py` has a `_load_family_outgroups()` function at line 36, but it is a separate script and not relevant to this file.)

## Test Note

**No live execution was possible in this environment.** The system's `python` on PATH is an old MGLTools 2.7 installation; no Python 3.11 environment with the required dependencies (Click, PyYAML, Biopython, ete3) is available. The fix was verified by:
1. Manual reading and parsing of the species.yaml schema
2. Counting entries under each section against the loader logic
3. Verifying line-by-line logic against the correct pattern in `root_tree.py`

Actual runtime testing must be performed in the conda carnivorenzyme environment on Hellbender HPC or a machine with Python 3.11 + dependencies installed.

## Verification Checklist

- [x] Original bug location: lines 37–49 of detect_convergence.py
- [x] Root cause identified: wrong YAML key ("species" vs "carnivorous_species"/"outgroup_species")
- [x] Fix applied correctly: iterates both sections, infers carnivorous status from membership
- [x] Return type unchanged: still returns {species_code: {carnivorous, carnivory_origin, name}}
- [x] Manual trace with real numbers: 18 carnivorous + 5 outgroup = 23 total ✓
- [x] No dead code found in this file
- [x] Pattern follows published reference in root_tree.py
- [x] No live execution environment available (noted explicitly)

## Verification

Independent re-check performed by a separate checker agent. No live Python execution was available in this environment either (confirmed: no working Python 3 with Click/PyYAML/Biopython/ete3 on PATH) — this verification is a manual code trace against the live files, not a test run.

### 1. Current source of `_load_species_meta()` (re-read fresh from disk)

Confirmed `workflow/scripts/detect_convergence.py` lines 37–61 now iterate `data.get("carnivorous_species", {}).items()` (setting `carnivorous=True` unconditionally, `carnivory_origin` from `entry.get("carnivory_origin", 0)`, `name=species_name`) and separately `data.get("outgroup_species", {}).items()` (setting `carnivorous=False` unconditionally, `carnivory_origin=0` hardcoded, `name=species_name`). The old `data.get("species", [])` pattern is gone from this function. Grep of the whole file for `data.get("species` / `"species"` key access confirms **no remaining occurrence** anywhere in `detect_convergence.py` — no dead code left over.

### 2. Downstream call sites of `species_meta` — all checked

Grepped every use of `species_meta` in the file:
- `_get_carnivorous_branch_lengths` (line 171): `species_meta.get(species, {}).get("carnivorous")` — safe `.get()`, matches new shape.
- `_detect_convergent_sites` (line 217): `known_species = set(species_meta.keys())` — keys are still species **codes**, unchanged contract.
- Line 248–251: `sm = species_meta.get(species, {})`; `sm.get("carnivorous")`; `sm.get("carnivory_origin", 0)` — all safe `.get()` calls, consistent with `{carnivorous, carnivory_origin, name}` shape.
- `main()` CLI block (line 457): `n_carni = sum(1 for v in species_meta.values() if v["carnivorous"])` — uses direct `v["carnivorous"]` indexing (not `.get`), but this is safe because **every** entry produced by the fixed loader (both carnivorous and outgroup branches) unconditionally sets the `"carnivorous"` key, so no `KeyError` risk.

No call site anywhere in the file still assumes the old (broken) shape or a top-level `species:` list. Confirmed clean.

### 3. `carnivorous` inferred from section membership, not an explicit field

Confirmed by reading both the code and the live YAML: the `carnivorous_species:` entries in `config/species.yaml` do **not** carry an explicit `carnivorous: true` field (e.g. `Nepenthes_gracilis` has no such key), yet the loader correctly assigns `True` purely because the entry came from the `carnivorous_species` section. Conversely, `outgroup_species:` entries (e.g. `Arabidopsis_thaliana`) *do* happen to carry an explicit `carnivorous: false` field in the YAML, but the loader ignores it and hardcodes `False` from section membership instead — which is the right bug-avoidance mechanism the audit doc claims (would still work even if that field were absent/inconsistent).

### 4. `carnivory_origin` pass-through

Confirmed the loader does not hardcode origin values for carnivorous species — it reads `entry.get("carnivory_origin", 0)` per-entry, so it will correctly reflect whatever origin numbers are in `config/species.yaml` even after those numbers are changed by a separate task. Outgroup species correctly get a hardcoded `0` (they have no `carnivory_origin` field in the YAML at all — confirmed by reading the outgroup section).

### 5. Independent recount against the live `config/species.yaml`

Manually enumerated every top-level key under each section by reading the full file (348 lines):

**`carnivorous_species:`** (18 entries, all with no explicit `carnivorous` field — status must come from section membership):
Nepenthes_gracilis (Ngra, origin 1), Nepenthes_alata (Nala, origin 1), Nepenthes_mirabilis (Nmir, origin 1), Drosera_capensis (Dcap, origin 2), Drosera_adelae (Dade, origin 2), Drosera_spatulata (Dspa, origin 2), Drosera_binata (Dbin, origin 2), Dionaea_muscipula (Dmus, origin 3), Aldrovanda_vesiculosa (Aves, origin 3), Cephalotus_follicularis (Cfol, origin 4), Sarracenia_purpurea (Spur, origin 5), Sarracenia_alata (Sala, origin 5), Darlingtonia_californica (Dcal, origin 5), Heliamphora_ciliata (Hcil, origin 5), Utricularia_gibba (Ugib, origin 6), Pinguicula_moranensis (Pmor, origin 6), Pinguicula_vulgaris (Pvul, origin 6), Byblis_filifolia (Bfil, origin 7).
**Count: 18.**

**`outgroup_species:`** (5 entries):
Arabidopsis_thaliana (Atha), Vitis_vinifera (Vvin), Solanum_lycopersicum (Slyo — note: code is genuinely `Slyo` not `Slyc` in the file, unusual but correctly transcribed), Hordeum_vulgare (Hvul), Secale_cereale (Scer).
**Count: 5.**

**Total: 23.**

**Result: matches the fixer's claimed "18 carnivorous + 5 outgroup = 23 total" exactly. No discrepancy found.** All individual codes and `carnivory_origin` values I independently transcribed also match the fixer's table entry-for-entry.

With this fix, `n_carni` (line 457) would report `18` and total species count `23` on any real run — not `0` as under the old bug.

### 6. Dead code check

Grepped the full file for any remaining reference to `.get("species"` / a top-level `species:` list pattern — none found. Also spot-checked `workflow/scripts/root_tree.py` (cited by the fixer as the reference pattern): it independently confirms `config.get("outgroup_species", {}).keys()` is indeed the established pattern elsewhere in this codebase, supporting the fixer's rationale (though root_tree.py only consumes `outgroup_species`, not `carnivorous_species`, for its own narrower purpose — not a discrepancy, just a narrower use of the same section-based pattern).

### Verdict: **CONFIRMED**

The fix correctly resolves the original bug. `_load_species_meta()` now reads `carnivorous_species` and `outgroup_species` instead of the nonexistent `species` key, infers `carnivorous` status purely from section membership (not a possibly-missing/inconsistent explicit field), passes through `carnivory_origin` per-entry without hardcoding, and every downstream consumer of `species_meta` in the file (`_get_carnivorous_branch_lengths`, `_detect_convergent_sites`, and the `main()` CLI logging) already expects and correctly handles the `{code: {carnivorous, carnivory_origin, name}}` shape. Independent recount of the live `config/species.yaml` reproduces the fixer's claimed 18 carnivorous + 5 outgroup = 23 total exactly, with no discrepancies in individual codes or origin values. No dead code or leftover references to the old broken key remain. This verification was performed as a manual trace only — no Python execution was possible in this environment — and that limitation is explicitly noted here, matching the fixer's own caveat.

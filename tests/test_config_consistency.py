"""Config-consistency tests for CarnivorEnzyme's enzyme_families.yaml / species.yaml.

Guards against the class of structural config bugs found in
audit/10_config_gaps_and_consistency_test.md, plus the blind spots closed in
audit/14_test_blindspots_and_doc_drift.md:
  (a) a species referenced as an accessions: key in enzyme_families.yaml that has no
      matching entry in species.yaml (silently drops that species from convergence
      detection — no error, just missing data)
  (b) a carnivorous_species entry with carnivory_origin unset/0 (silently excluded from
      every cross-lineage convergence count, since detect_convergence.py skips origin==0)
  (c) a family key defined in both tier1: and methods_benchmark: (ambiguous which
      section governs it; Snakefile's FAMILIES list only reads tier1:, so a family
      present in both would silently behave as tier1 only, masking a methods_benchmark
      intent)
  (d) the species-name-vs-code keying bug in detect_convergence.py's _seq_id_to_species():
      tree leaves / FASTA headers are named "{Species_name}|{accession}" by
      fetch_sequences.py, so species.yaml lookups must key on the full species name
      (not just the short `code`) for the mapping to succeed at all.
  (e) a species listed in BOTH carnivorous_species: and outgroup_species: — since
      detect_convergence.py's _load_species_meta() loads carnivorous_species first then
      outgroup_species second, the second silently overwrites the first in its dict,
      reclassifying a carnivorous species as non-carnivorous (carnivory_origin=0). The
      union-based `_known_species_keys()` helper used elsewhere in this file cannot catch
      this (a key defined twice is indistinguishable from a key defined once in a set
      union), so this needs its own dedicated disjointness check.
  (f) duplicate `code:` values across carnivorous_species: and outgroup_species: combined
      — `_load_species_meta()` also keys its dict by `code`, so two species sharing a
      `code` would have the second silently overwrite the first the same way as (e).
  (g) a tier1: family whose non-TODO accessions span fewer than
      config.yaml's convergence.min_lineages distinct carnivory_origin values — this is
      exactly the invariant that motivated moving nepenthesins/neprosins and
      chitinases_gh19_class_i out of tier1 in earlier remediation (see
      audit/03_merops_restructure_and_neprosin_rescope.md,
      audit/10_config_gaps_and_consistency_test.md); it can silently regress if someone
      edits accessions later without noticing.

Real PyYAML parsing of the actual config files — no mocking of config content.
"""

import importlib.util
import sys
import types
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).parent.parent
ENZYME_FAMILIES_YAML = PROJECT_ROOT / "config" / "enzyme_families.yaml"
SPECIES_YAML = PROJECT_ROOT / "config" / "species.yaml"
CONFIG_YAML = PROJECT_ROOT / "config" / "config.yaml"
DETECT_CONVERGENCE_PY = PROJECT_ROOT / "workflow" / "scripts" / "detect_convergence.py"


# ---------------------------------------------------------------------------
# Shared loading helpers
# ---------------------------------------------------------------------------

def _load_enzyme_families() -> dict:
    with open(ENZYME_FAMILIES_YAML, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _load_species() -> dict:
    with open(SPECIES_YAML, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _load_config() -> dict:
    with open(CONFIG_YAML, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _accession_species_keys(enzyme_families: dict) -> set[str]:
    """Every species key used under any family's accessions: mapping (tier1 + methods_benchmark)."""
    used = set()
    for section in ("tier1", "methods_benchmark"):
        for _family_name, family_data in enzyme_families.get(section, {}).items():
            accessions = family_data.get("accessions", {}) or {}
            used.update(accessions.keys())
    return used


def _known_species_keys(species: dict) -> set[str]:
    """Every species name in species.yaml's carnivorous_species: and outgroup_species: sections."""
    return set(species.get("carnivorous_species", {}).keys()) | set(
        species.get("outgroup_species", {}).keys()
    )


def _import_detect_convergence():
    """Import workflow/scripts/detect_convergence.py as a module, by file path.

    It isn't a package (no __init__.py) and isn't installed, so it can't be imported
    by dotted name. It also does `import ete3` at module level; ete3 is a heavy
    optional dependency that may not be installed in every environment running this
    test suite (it wasn't in the environment this test was authored/verified in). We
    only need the pure-Python `_seq_id_to_species` / `_load_species_meta` helpers,
    which don't touch ete3 at all, so stub it out if the real package is missing
    rather than skip the test.
    """
    if "ete3" not in sys.modules:
        try:
            import ete3  # noqa: F401
        except ImportError:
            stub = types.ModuleType("ete3")

            class _StubTree:  # pragma: no cover - never instantiated by these tests
                def __init__(self, *args, **kwargs):
                    raise NotImplementedError(
                        "ete3 is stubbed out for test_config_consistency.py; "
                        "install the real ete3 package to exercise tree-parsing code paths."
                    )

            stub.Tree = _StubTree
            sys.modules["ete3"] = stub

    spec = importlib.util.spec_from_file_location(
        "detect_convergence_under_test", DETECT_CONVERGENCE_PY
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# (a) every accession species key exists in species.yaml
# ---------------------------------------------------------------------------

def test_all_accession_species_exist_in_species_yaml():
    enzyme_families = _load_enzyme_families()
    species = _load_species()

    used = _accession_species_keys(enzyme_families)
    known = _known_species_keys(species)

    missing = sorted(used - known)
    assert not missing, (
        "Species referenced under a family's accessions: mapping in "
        f"{ENZYME_FAMILIES_YAML.name} but absent from both carnivorous_species: and "
        f"outgroup_species: in {SPECIES_YAML.name}: {missing}"
    )


# ---------------------------------------------------------------------------
# (b) every carnivorous_species entry has a nonzero carnivory_origin
# ---------------------------------------------------------------------------

def test_carnivorous_species_have_nonzero_carnivory_origin():
    species = _load_species()
    carnivorous = species.get("carnivorous_species", {})
    assert carnivorous, "species.yaml carnivorous_species: is empty or missing"

    bad = {
        name: entry.get("carnivory_origin")
        for name, entry in carnivorous.items()
        if not entry.get("carnivory_origin")
    }
    assert not bad, (
        "carnivorous_species entries with a missing/zero carnivory_origin (these are "
        f"silently excluded from convergence detection): {bad}"
    )


# ---------------------------------------------------------------------------
# (c) no family key defined in both tier1: and methods_benchmark:
# ---------------------------------------------------------------------------

def test_tier1_and_methods_benchmark_families_disjoint():
    enzyme_families = _load_enzyme_families()
    tier1_keys = set(enzyme_families.get("tier1", {}).keys())
    methods_benchmark_keys = set(enzyme_families.get("methods_benchmark", {}).keys())

    overlap = tier1_keys & methods_benchmark_keys
    assert overlap == set(), (
        f"Family key(s) defined under BOTH tier1: and methods_benchmark: {sorted(overlap)} "
        "— Snakefile derives FAMILIES only from tier1:, so this would silently drop the "
        "methods_benchmark definition from convergence targets while leaving it ambiguous "
        "which block's metadata (accessions, notes) is authoritative."
    )


# ---------------------------------------------------------------------------
# (d) _seq_id_to_species round-trip on every species.yaml name
# ---------------------------------------------------------------------------

def test_seq_id_to_species_roundtrip_on_species_names():
    """fetch_sequences.py writes FASTA/tree-leaf IDs as '{Species_name}|{accession}'
    (record.id = f"{species}|{acc}", species = the full name key from enzyme_families.yaml
    accessions:, which is the same string used as the top-level key in species.yaml).
    detect_convergence.py's _seq_id_to_species() must map such an ID back to that same
    species name when given the ACTUAL known_species set produced by
    _load_species_meta() — this is exactly the species-name-vs-code keying bug this test
    guards against.

    NOTE (audit/14_test_blindspots_and_doc_drift.md): a prior version of this test built
    its `known_species` set independently from species.yaml directly (via
    `_known_species_keys()`), rather than from `_load_species_meta()`'s actual output.
    That meant it could NOT catch a regression in the loader itself (e.g. if
    `_load_species_meta` reverted to keying only by short `code`, dropping the full
    species-name keys) — it was only exercising `_seq_id_to_species()` in isolation, not
    the loader it claimed to guard per its own docstring. Fixed to derive `known_species`
    from `set(module._load_species_meta(SPECIES_YAML).keys())` so it actually exercises
    the loader. (`test_load_species_meta_keys_include_full_species_names` below already
    separately guards the loader's key set directly, so this was not a coverage gap for
    the loader — just a docstring that overclaimed what this specific test exercised.)
    """
    module = _import_detect_convergence()
    species = _load_species()
    species_names_to_roundtrip = _known_species_keys(species)
    assert species_names_to_roundtrip, "no species names found in species.yaml to round-trip"

    known_species = set(module._load_species_meta(SPECIES_YAML).keys())

    failures = []
    for species_name in sorted(species_names_to_roundtrip):
        fake_seq_id = f"{species_name}|FAKEACC.1"
        mapped = module._seq_id_to_species(fake_seq_id, known_species)
        if mapped != species_name:
            failures.append((fake_seq_id, mapped))

    assert not failures, f"_seq_id_to_species round-trip failures (seq_id -> mapped): {failures}"


def test_load_species_meta_keys_include_full_species_names():
    """_load_species_meta's returned dict must be keyed by full species names (not only
    short codes), since that's what tree leaves / FASTA headers actually use. This is a
    direct regression guard for the species-name-vs-code keying bug fixed in
    detect_convergence.py._load_species_meta (registers each species under both its
    species.yaml key and its `code`).
    """
    module = _import_detect_convergence()
    meta = module._load_species_meta(SPECIES_YAML)

    species = _load_species()
    known_species_names = _known_species_keys(species)

    missing = sorted(known_species_names - set(meta.keys()))
    assert not missing, (
        f"_load_species_meta({SPECIES_YAML.name}) is missing full-species-name keys for: "
        f"{missing} — only short `code` keys were registered for these, which breaks "
        "matching against '{Species_name}|{accession}'-style tree leaf names."
    )


# ---------------------------------------------------------------------------
# (e) no species listed in BOTH carnivorous_species: and outgroup_species:
# ---------------------------------------------------------------------------

def test_carnivorous_and_outgroup_species_disjoint():
    """No species name may appear as a key in both carnivorous_species: and
    outgroup_species:.

    _load_species_meta() in detect_convergence.py loads carnivorous_species first, then
    outgroup_species second, with each pass unconditionally overwriting
    `meta[species_name]`. A species listed in both sections would therefore be silently
    reclassified as non-carnivorous (carnivorous=False, carnivory_origin=0) by whichever
    section is loaded last — dropping it from every cross-lineage convergence count with
    no error. The union-based `_known_species_keys()` helper used by other tests in this
    file cannot detect this: a key present in both sets of a Python set union is
    indistinguishable from a key present in only one.
    """
    species = _load_species()
    carnivorous_keys = set(species.get("carnivorous_species", {}).keys())
    outgroup_keys = set(species.get("outgroup_species", {}).keys())

    overlap = carnivorous_keys & outgroup_keys
    assert overlap == set(), (
        f"Species listed in BOTH carnivorous_species: and outgroup_species: in "
        f"{SPECIES_YAML.name}: {sorted(overlap)} — detect_convergence.py's "
        "_load_species_meta() loads carnivorous_species first then outgroup_species "
        "second, so the second entry silently overwrites the first, reclassifying this "
        "species as non-carnivorous with carnivory_origin=0."
    )


# ---------------------------------------------------------------------------
# (f) no duplicate `code:` values across carnivorous_species: + outgroup_species:
# ---------------------------------------------------------------------------

def test_species_codes_are_unique():
    """No two species (across carnivorous_species: and outgroup_species: combined) may
    share the same `code:` value.

    _load_species_meta() also registers each species under `meta[code]`, in the same
    load order as the name keys (carnivorous_species first, then outgroup_species). Two
    species sharing a `code` would have the second silently overwrite the first's `code`
    entry in that dict — the same class of silent-data-loss bug as (e), but keyed on
    `code` instead of the species name.
    """
    species = _load_species()
    seen: dict[str, str] = {}
    dupes: dict[str, list[str]] = {}
    for section in ("carnivorous_species", "outgroup_species"):
        for name, entry in species.get(section, {}).items():
            code = entry.get("code")
            if not code:
                continue
            if code in seen:
                dupes.setdefault(code, [seen[code]]).append(name)
            else:
                seen[code] = name

    assert not dupes, (
        f"Duplicate `code:` values in {SPECIES_YAML.name} (code -> species names): "
        f"{dupes} — _load_species_meta() keys its dict by `code` as well as by species "
        "name, so the second species sharing a code silently overwrites the first."
    )


# ---------------------------------------------------------------------------
# (g) every tier1: family spans >= convergence.min_lineages carnivory_origin values
#     among its non-TODO accessions
# ---------------------------------------------------------------------------

def _is_todo_only(accession_entry) -> bool:
    """True if a family's accessions: entry for one species has no real accession.

    Entries are lists of accession strings; a placeholder entry is `[TODO]` (a bare
    'TODO' sentinel, case-insensitive) with no real accession alongside it.
    """
    if not isinstance(accession_entry, list):
        return True
    real = [a for a in accession_entry if isinstance(a, str) and a.strip().upper() != "TODO"]
    return not real


def test_tier1_families_span_min_lineages_carnivory_origins():
    """Every tier1: family must have non-TODO accessions from at least
    config.yaml's convergence.min_lineages distinct carnivory_origin values.

    A tier1 family confined to a single carnivory_origin can mathematically never
    produce a cross-lineage convergence result under detect_convergence.py's
    min_lineages gate (see _detect_convergent_sites: `if n_origins < min_lineages:
    continue`) — every run silently produces a header-only empty output. This is exactly
    why chitinases_gh19_class_i, nepenthesins, and neprosins were moved out of tier1:
    into methods_benchmark: in earlier remediation (audit/03 and audit/10). This test
    guards against that regressing if someone edits accessions later (e.g. TODO-ing out
    the last non-TODO species from a second origin) without noticing the family no
    longer qualifies for tier1:.
    """
    enzyme_families = _load_enzyme_families()
    species = _load_species()
    config = _load_config()

    min_lineages = config["convergence"]["min_lineages"]

    origin_by_species = {
        name: entry.get("carnivory_origin")
        for name, entry in species.get("carnivorous_species", {}).items()
    }

    failures = {}
    for family, family_data in enzyme_families.get("tier1", {}).items():
        accessions = family_data.get("accessions", {}) or {}
        origins = set()
        for species_name, accession_entry in accessions.items():
            if species_name not in origin_by_species:
                continue  # not a carnivorous species (e.g. an outgroup) — doesn't count
            if _is_todo_only(accession_entry):
                continue
            origins.add(origin_by_species[species_name])

        if len(origins) < min_lineages:
            failures[family] = sorted(origins)

    assert not failures, (
        f"tier1: family/families with non-TODO accessions spanning fewer than "
        f"min_lineages={min_lineages} distinct carnivory_origin values (family -> "
        f"origins found): {failures} — such a family can never produce a cross-lineage "
        "convergence result and should be moved to methods_benchmark: (see "
        "audit/03_merops_restructure_and_neprosin_rescope.md and "
        "audit/10_config_gaps_and_consistency_test.md for the precedent)."
    )

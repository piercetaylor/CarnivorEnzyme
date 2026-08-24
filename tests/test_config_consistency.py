"""Config-consistency tests for CarnivorEnzyme's enzyme_families.yaml / species.yaml.

Guards against the class of structural config bugs found in
audit/10_config_gaps_and_consistency_test.md:
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
    species name when given a known_species set containing species names (as
    _load_species_meta's returned dict keys do, alongside short codes) — this is exactly
    the species-name-vs-code keying bug this test guards against.
    """
    module = _import_detect_convergence()
    species = _load_species()
    known_species_names = _known_species_keys(species)
    assert known_species_names, "no species names found in species.yaml to round-trip"

    failures = []
    for species_name in sorted(known_species_names):
        fake_seq_id = f"{species_name}|FAKEACC.1"
        mapped = module._seq_id_to_species(fake_seq_id, known_species_names)
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

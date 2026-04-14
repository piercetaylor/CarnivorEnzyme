"""CarnivorEnzyme master Snakemake workflow.

Execution order:
  Phase 1  — Sequence retrieval          (retrieve.smk)
  Phase 2  — Alignment + Phylogeny       (alignment.smk, phylogeny.smk)
  Phase 2+ — Orthology (optional)        (orthology.smk)
  Phase 3  — Convergence detection       (convergence.smk)
  Phase 4  — Structure prediction (HPC)  (predict_structure.smk, structural_align.smk)
  Phase 5A — FoldX stability             (foldx.smk)
  Phase 5B — EVcouplings/EVmutation      (evmutation.smk)
  Phase 5C — FoldX × EVmutation compare  (integrate.smk)
  Phase 6  — Molecular docking           (docking.smk)
  Phase 7  — Electrostatics              (electrostatics.smk)
  Phase 8  — Expression (optional)       (expression.smk)
  Phase 9  — Atlas + Figures             (integrate.smk)

Usage:
  # Run full pipeline (all families)
  snakemake --use-conda --cores 8

  # Run Phase 1 only (sequence download)
  snakemake --use-conda phase1 --cores 1

  # Run Phase 1+2 for a specific family
  snakemake --use-conda --cores 8 \\
      results/alignments/nepenthesins.trimmed.fa \\
      results/phylogenies/nepenthesins.treefile

  # Dry-run to check DAG
  snakemake --use-conda -n
"""

import yaml
from pathlib import Path

configfile: "config/config.yaml"

# ── Derive FAMILIES list from enzyme_families.yaml at parse time ──────────────
def _load_tier1_families(path: str) -> list[str]:
    with open(path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    return list(data.get("tier1", {}).keys())

FAMILIES = _load_tier1_families("config/enzyme_families.yaml")


# ── Include all rule modules ───────────────────────────────────────────────────
include: "workflow/rules/retrieve.smk"
include: "workflow/rules/orthology.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/phylogeny.smk"
include: "workflow/rules/convergence.smk"
include: "workflow/rules/predict_structure.smk"
include: "workflow/rules/structural_align.smk"
include: "workflow/rules/foldx.smk"
include: "workflow/rules/evmutation.smk"
include: "workflow/rules/docking.smk"
include: "workflow/rules/electrostatics.smk"
include: "workflow/rules/expression.smk"
include: "workflow/rules/integrate.smk"


# ── Phase-level convenience targets ───────────────────────────────────────────

rule phase1:
    """Phase 1 complete: all sequences downloaded."""
    input:
        "results/sequences/",


rule phase2:
    """Phase 2 complete: trimmed alignments and ML trees for all families."""
    input:
        trimmed=expand(
            "results/alignments/{family}.trimmed.fa",
            family=FAMILIES,
        ),
        trees=expand(
            "results/phylogenies/{family}.treefile",
            family=FAMILIES,
        ),
        rooted=expand(
            "results/phylogenies/{family}.rooted.treefile",
            family=FAMILIES,
        ),


# ── Final target ──────────────────────────────────────────────────────────────

rule all:
    """Full pipeline: structural atlas SQLite + all manuscript figures."""
    input:
        "results/atlas/atlas.sqlite",
        expand("results/atlas/figures/fig{n}.pdf", n=range(1, 9)),

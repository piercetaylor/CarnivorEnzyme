"""CarnivorEnzyme master Snakemake workflow.

Execution order:
  Phase 1   — Sequence retrieval           (retrieve.smk)
  Phase 2   — Alignment + Phylogeny        (alignment.smk, phylogeny.smk)
  Phase 2+  — Orthology (optional)         (orthology.smk)
  Phase 3   — Convergence detection        (convergence.smk)
  Phase 3A  — ESM-2 fitness triage         (esm2.smk)           ← NEW
  Phase 4   — Structure prediction (HPC)   (predict_structure.smk, structural_align.smk)
  Phase 4A  — Ancestral reconstruction     (ancestral_structure.smk) ← NEW
  Phase 5A  — FoldX stability              (foldx.smk)
  Phase 5B  — EVcouplings/EVmutation       (evmutation.smk)
  Phase 5C  — FoldX × EVmutation compare   (integrate.smk)
  Phase 5D  — Alchemical FEP binding ΔΔG   (fep.smk)            ← NEW
  Phase 5E  — Constant-pH MD (CpHMD)       (cphmd.smk)          ← NEW
  Phase 5F  — QM/MM (gated, ≤2 targets)   (qmmm.smk)           ← NEW (disabled by default)
  Phase 6   — Molecular docking            (docking.smk)
  Phase 7   — Electrostatics               (electrostatics.smk)
  Phase 8   — Expression (optional)        (expression.smk)
  Phase 9   — Atlas + Figures              (integrate.smk)

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
include: "workflow/rules/esm2.smk"                  # Phase 3A: ESM-2 triage
include: "workflow/rules/predict_structure.smk"
include: "workflow/rules/ancestral_structure.smk"    # Phase 4A: ancestral reconstruction
include: "workflow/rules/structural_align.smk"
include: "workflow/rules/foldx.smk"
include: "workflow/rules/evmutation.smk"
include: "workflow/rules/fep.smk"                   # Phase 5D: alchemical FEP
include: "workflow/rules/cphmd.smk"                 # Phase 5E: constant-pH MD
include: "workflow/rules/qmmm.smk"                  # Phase 5F: QM/MM (disabled by default)
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


rule phase3a:
    """Phase 3A complete: ESM-2 LLR scores for convergent sites in all families."""
    input:
        expand("results/esm2/{family}.esm2_scores.tsv", family=FAMILIES),


rule phase4a:
    """Phase 4A complete: ancestral sequences extracted + structures predicted."""
    input:
        ancestors=expand(
            "results/ancestral/{family}.mrca_ancestor.fa",
            family=FAMILIES,
        ),
        ancestor_pdbs=expand(
            "results/ancestral/{family}.mrca_ancestor.pdb",
            family=FAMILIES,
        ),
        comparisons=expand(
            "results/ancestral/{family}.structural_comparison.tsv",
            family=FAMILIES,
        ),


rule phase5d:
    """Phase 5D complete: FEP ΔΔG_bind results for all families."""
    input:
        expand("results/fep/{family}.fep_results.tsv", family=FAMILIES),


# ── Final target ──────────────────────────────────────────────────────────────

rule all:
    """Full pipeline: structural atlas SQLite + all manuscript figures.

    Figures:
      fig1  — Per-family phylogenies with convergent sites on branches
      fig2  — Structural superposition per family (convergent sites highlighted)
      fig3  — FoldX ΔΔG heatmap (positions × lineages)
      fig4  — EVmutation ΔΔE heatmap (same positions × lineages)
      fig5  — FoldX × EVmutation quadrant scatter plots per family  [core novel figure]
      fig6  — ESM-2 LLR vs FoldX ΔΔG correlation
      fig7  — Docking affinity comparison across orthologs
      fig8  — Electrostatic surfaces at pH 2.5 vs 5.0
      fig9  — FEP ΔΔG_bind for priority convergent sites
      fig10 — Ancestral vs modern active-site geometry (RMSD_convergent)
    """
    input:
        "results/atlas/atlas.sqlite",
        expand("results/atlas/figures/fig{n}.pdf", n=range(1, 11)),

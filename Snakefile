"""CarnivorEnzyme master Snakemake workflow."""

configfile: "config/config.yaml"

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


FAMILIES = list(config.get("families_to_run", ["neprosins"]))


rule all:
    input:
        "results/atlas/atlas.sqlite",
        expand("results/atlas/figures/fig{n}.pdf", n=range(1, 9)),

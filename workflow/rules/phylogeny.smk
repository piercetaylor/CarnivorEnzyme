"""phylogeny.smk — Maximum-likelihood phylogeny + ancestral reconstruction (Phase 2).

TWO IQ-TREE PASSES, not one. See audit/11_ancestral_reconstruction_architecture_fix.md.

  1. infer_tree               ML topology search, ModelFinder Plus, UFBoot2 + SH-aLRT.
                              NO -ancestral here.
                              -> results/phylogenies/{family}.treefile   (support labels)
                              -> results/phylogenies/{family}.iqtree     (ModelFinder verdict)

  2. root_tree                Outgroup-root that topology with root_tree.py.
                              -> {family}.rooted.treefile        (for human inspection)
                              -> {family}.rooted.plain.treefile  (labels stripped, for -te)
                              -> {family}.outgroups.txt          (tip labels, for -o)

  3. ancestral_reconstruction Re-run IQ-TREE with the rooted topology FIXED (-te) and
                              --ancestral, reusing pass 1's ModelFinder verdict.
                              -> {family}.asr.treefile   (NodeN internal labels)
                              -> {family}.asr.state      (Node column == those labels)

Why the second pass exists
--------------------------
IQ-TREE's `.state` file refers to the internal-node labels of the treefile that the SAME
invocation wrote — its own header says "Ancestral state reconstruction for all nodes in
<prefix>.treefile". Running `-ancestral` during pass 1 therefore produces states keyed to
the UNROOTED pass-1 tree, while everything downstream consumed the separately rooted tree.
Rerooting moves those labels onto different clades and can add nodes with no `.state`
entry, so ancestral states were being read off the wrong nodes. Reconstructing AFTER
rooting, and consuming IQ-TREE's own output tree, removes the correspondence problem
instead of patching around it.

Verified against IQ-TREE 3.1.3 (Windows x86-64 release binary), 2026-08-24:
  * `-te` fixes the topology (0 search iterations) and is accepted for a rooted Newick;
    IQ-TREE logs "rooted tree" and then silently un-roots it, which is why the rooted
    tree's own node structure must never be assumed to match `.state`.
  * `-B`/`-bb` is REJECTED alongside `-te`: "Ultrafast bootstrap does not work with
    -fast, -te or -n option". Support values belong to pass 1 only.
  * Feeding a label-stripped tree makes IQ-TREE write clean `Node1..NodeK` labels into
    `<prefix>.treefile`, matching `<prefix>.state` exactly.
  * Marginal ancestral states are identical (max |Δp| = 1e-5, i.e. output rounding)
    between the pass-1 unrooted run and the `-te` rooted run at every shared vertex,
    as expected for a reversible model.

NOTE on binary name:
  IQ-TREE 3 installs as 'iqtree3' in the conda package (bioinfo.yaml pins iqtree>=3.0.0).
  If the cluster module exposes it as 'iqtree', override config["phylogeny"]["binary"].
"""


rule infer_tree:
    """Infer the ML phylogeny with IQ-TREE 3 (topology search pass).

    Key flags:
      -m MFP        ModelFinder Plus: automatically selects best-fit model
      -bb 1000      UFBoot2 ultrafast bootstrap (1000 replicates)
      -alrt 1000    SH-aLRT branch test (reported alongside UFBoot2)
      -seed 42      Reproducible random seed
      -redo         Overwrite previous run in same prefix directory
      -pre          Output prefix (all IQ-TREE files share this stem)

    -ancestral is deliberately NOT used here; ancestral states are reconstructed on the
    rooted topology by the ancestral_reconstruction rule below.

    UFBoot2 interpretation: values >=95 indicate strong support (NOT comparable
    to standard bootstrap >=70 cutoff). Use alrt >=80 AND ufboot >=95 together.
    """
    input:
        alignment="results/alignments/{family}.trimmed.fa",
    output:
        treefile="results/phylogenies/{family}.treefile",
        iqlog="results/phylogenies/{family}.iqtree",
    params:
        bootstrap=config["phylogeny"]["bootstrap"],
        model=config["phylogeny"]["model"],
        binary=config["phylogeny"].get("binary", "iqtree3"),
        prefix="results/phylogenies/{family}",
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/phylogeny/{family}.log",
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240,
    shell:
        """
        mkdir -p results/phylogenies/
        {params.binary} -s {input.alignment} -m {params.model} \
            -bb {params.bootstrap} -alrt 1000 -nt {threads} -seed 42 \
            -redo -pre {params.prefix} 2> {log}
        grep "Best-fit model:" {output.iqlog} >> {log} 2>/dev/null || true
        """


rule root_tree:
    """Root the IQ-TREE topology on the outgroup species listed in species.yaml.

    Emits three files: the rooted tree (human/figure use), the same tree with internal
    labels stripped (IQ-TREE `-te` input — pass-1 support labels would otherwise be
    copied verbatim into the `.state` Node column while describing the wrong clades),
    and the matched outgroup tip labels (IQ-TREE `-o` input).
    """
    input:
        treefile="results/phylogenies/{family}.treefile",
        species_config="config/species.yaml",
        families_config="config/enzyme_families.yaml",
    output:
        rooted="results/phylogenies/{family}.rooted.treefile",
        plain="results/phylogenies/{family}.rooted.plain.treefile",
        outgroups="results/phylogenies/{family}.outgroups.txt",
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/phylogeny/{family}.root.log",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10,
    shell:
        """
        python workflow/scripts/root_tree.py \
            --treefile {input.treefile} \
            --species-config {input.species_config} \
            --families-config {input.families_config} \
            --family {wildcards.family} \
            --output {output.rooted} \
            --output-plain {output.plain} \
            --output-outgroups {output.outgroups} \
            --verbose \
            2> {log}
        """


rule ancestral_reconstruction:
    """Marginal ancestral state reconstruction on the FIXED, already-rooted topology.

    run_ancestral.py reuses pass 1's ModelFinder verdict (parsed from {family}.iqtree),
    runs `iqtree3 -te <rooted, label-stripped> -o <outgroups> --ancestral`, and then
    verifies that the emitted .asr.treefile and .asr.state describe exactly the same node
    set and that the topology was not altered. Any mismatch raises.
    """
    input:
        alignment="results/alignments/{family}.trimmed.fa",
        plain="results/phylogenies/{family}.rooted.plain.treefile",
        outgroups="results/phylogenies/{family}.outgroups.txt",
        report="results/phylogenies/{family}.iqtree",
    output:
        asr_tree="results/phylogenies/{family}.asr.treefile",
        state="results/phylogenies/{family}.asr.state",
        asr_report="results/phylogenies/{family}.asr.iqtree",
    params:
        binary=config["phylogeny"].get("binary", "iqtree3"),
        fallback_model=config["phylogeny"].get("fallback_model", "LG+G4"),
        prefix="results/phylogenies/{family}.asr",
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/phylogeny/{family}.ancestral.log",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        """
        python workflow/scripts/run_ancestral.py \
            --alignment {input.alignment} \
            --tree {input.plain} \
            --report {input.report} \
            --outgroups {input.outgroups} \
            --prefix {params.prefix} \
            --iqtree-binary {params.binary} \
            --fallback-model {params.fallback_model} \
            --threads {threads} \
            --seed 42 \
            --verbose \
            2> {log}
        """

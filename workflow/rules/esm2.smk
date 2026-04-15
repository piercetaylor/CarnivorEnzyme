"""esm2.smk — Phase 3A: ESM-2 masked marginal scoring of convergent substitutions.

ESM-2 (Lin et al. 2023, Science 379:1123) provides a zero-shot fitness estimate
purely from sequence, independent of FoldX (structure-based) and EVmutation
(coevolution-based). The log-likelihood ratio (LLR) for each convergent substitution
acts as a triage filter:

  LLR > 0  → language model favors derived AA in this sequence context
  LLR ≤ 0  → model does not support the derived AA; convergence may be neutral or false

This rule runs BEFORE FoldX (Phase 5A) so high-LLR sites can be prioritized for
expensive downstream analyses (CpHMD, FEP).

Input:
  results/convergence/{family}.convergent_sites.tsv  (from convergence.smk)
  results/sequences/{family}/all_sequences.fa         (concatenated per-family FASTA)

Output:
  results/esm2/{family}.esm2_scores.tsv

Citation:
  Lin Z, Akin H, Rao R, et al. (2023) Evolutionary-scale prediction of atomic-level
  protein structure with a language model. Science 379:1123–1130.
  Meier J, Rao R, Verkuil R, et al. (2021) Language models enable zero-shot prediction
  of the effects of mutations on protein function. NeurIPS 34.
"""


rule score_esm2:
    """Score convergent substitutions with ESM-2 masked marginal log-likelihood.

    Runs once per enzyme family. Uses GPU if available (config.esm2.device).
    Model is downloaded on first run to ~/.cache/torch/hub/ (~1.5 GB for 650M model).
    """
    input:
        convergent="results/convergence/{family}.convergent_sites.tsv",
        sequences="results/sequences/{family}/all_sequences.fa",
    output:
        tsv="results/esm2/{family}.esm2_scores.tsv",
    params:
        model=config["esm2"]["model"],
        batch_size=config["esm2"]["batch_size"],
        device=config["esm2"]["device"],
    conda:
        "../envs/esm.yaml"
    log:
        "logs/esm2/{family}.log",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
        # GPU resource: Snakemake SLURM profile maps this to --gres=gpu:1
        gpu=1,
    shell:
        """
        python workflow/scripts/score_esm2.py \
            --convergent {input.convergent} \
            --sequences  {input.sequences} \
            --output     {output.tsv} \
            --model-name {params.model} \
            --batch-size {params.batch_size} \
            --device     {params.device} \
            --verbose \
            2> {log}
        """

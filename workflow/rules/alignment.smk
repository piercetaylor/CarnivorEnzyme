"""alignment.smk — Multiple sequence alignment and trimming (Phase 2).

Two rules per family:
  1. align_family   — MAFFT L-INS-i (--localpair --maxiterate 1000)
  2. trim_alignment — trimAl -automated1 (removes columns >50% gaps)

Input:  results/sequences/{family}/combined.fa  (from retrieve.smk)
Output: results/alignments/{family}.trimmed.fa  (used by phylogeny.smk)

Tool justification (CLAUDE.md §1):
  MAFFT L-INS-i: most accurate mode for <200 sequences; iterative local pairwise.
  trimAl -automated1: automated gap threshold preserving informative positions.
  Katoh & Standley 2013 (MBE 30:772); Capella-Gutiérrez et al. 2009 (Bioinf 25:1972).
"""


rule align_family:
    """Run MAFFT L-INS-i on the combined per-family FASTA.

    --localpair + --maxiterate 1000: most accurate alignment mode.
    Appropriate because each family has <100 sequences.
    """
    input:
        fasta="results/sequences/{family}/combined.fa",
    output:
        aligned="results/alignments/{family}.aligned.fa",
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/alignment/{family}.mafft.log",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        """
        mafft \
            --localpair \
            --maxiterate 1000 \
            --thread {threads} \
            --reorder \
            {input.fasta} \
            > {output.aligned} \
            2> {log}

        echo "Aligned $(grep -c '^>' {output.aligned}) sequences" >> {log}
        """


rule trim_alignment:
    """Trim gapped alignment columns with trimAl -automated1.

    -automated1: automatic heuristic that selects the gap threshold giving
    the highest average identity among retained columns. Suitable for
    downstream phylogenetics (does not over-trim).

    Also outputs a summary of retained / removed columns to the log.
    """
    input:
        aligned="results/alignments/{family}.aligned.fa",
    output:
        trimmed="results/alignments/{family}.trimmed.fa",
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/alignment/{family}.trimal.log",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15,
    shell:
        """
        trimal \
            -in {input.aligned} \
            -out {output.trimmed} \
            -automated1 \
            -htmlout results/alignments/{wildcards.family}.trimal.html \
            2> {log}

        echo "Trimmed alignment: $(grep -c '^>' {output.trimmed}) seqs" >> {log}

        # Count alignment columns before and after trimming
        python3 -c "
from Bio import SeqIO
raw = next(SeqIO.parse('{input.aligned}', 'fasta'))
trim = next(SeqIO.parse('{output.trimmed}', 'fasta'))
print(f'Columns before trimming: {{len(raw.seq)}}', flush=True)
print(f'Columns after  trimming: {{len(trim.seq)}}', flush=True)
print(f'Removed: {{len(raw.seq) - len(trim.seq)}} columns', flush=True)
" >> {log} 2>&1
        """

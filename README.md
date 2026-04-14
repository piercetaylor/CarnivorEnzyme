# CarnivorEnzyme
<<<<<<< HEAD

Structural atlas of convergently evolved digestive enzymes in carnivorous plants.

## What This Is

A reproducible Snakemake pipeline that predicts 3D structures of digestive enzymes from all independently-evolved carnivorous plant lineages, maps convergent amino acid substitutions onto those structures, quantifies thermodynamic and evolutionary consequences, and predicts substrate specificity shifts via molecular docking.

## Quick Start

```bash
# 1. Clone
git clone https://github.com/YOUR_USERNAME/CarnivorEnzyme.git
cd CarnivorEnzyme

# 2. Create conda environment
conda env create -f environment.yml
conda activate carnivorenzyme

# 3. Verify external tools
foldx --help          # Requires academic license from foldxsuite.crg.eu
which TMalign         # Compile from zhanggroup.org/TM-align

# 4. Configure
cp config/config.yaml.example config/config.yaml
# Edit config.yaml with your NCBI email, AF3 path, etc.

# 5. Dry run
snakemake -n

# 6. Run (local)
snakemake --cores 8 --use-conda

# 7. Run (Hellbender SLURM)
snakemake --profile config/slurm/
```

## Project Guide

See [CLAUDE.md](CLAUDE.md) for the complete project specification, method justifications, build order, and test gates. This file is also read by Claude Code for automated development.

## Citation

If you use this pipeline or database, please cite:

> Taylor P. (2026) CarnivorEnzyme: A Pan-Carnivorous Structural Atlas of Convergently Evolved Digestive Enzymes with Functional Prediction. *bioRxiv* [preprint].

## License

MIT
=======
Predicting 3D structures of digestive enzymes from independently-evolved carnivorous plant lineages
>>>>>>> 277cae1871acb45a5cc5d19b90df6d44aeb62242

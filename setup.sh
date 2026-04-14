#!/usr/bin/env bash
# setup.sh — Initialize the CarnivorEnzyme repository.
# Usage: bash setup.sh
set -euo pipefail

echo "=== CarnivorEnzyme: Setting up repository ==="

# ── Directories ──────────────────────────────────────────────────────────────
dirs=(
  config config/slurm
  workflow/rules workflow/scripts workflow/envs
  resources/accessions resources/structures/experimental resources/substrates
  results
  webapp/pages
  manuscript/figures manuscript/tables
  tests/data
)
for d in "${dirs[@]}"; do mkdir -p "$d"; done
echo "  ✓ directories created"

# ── .gitignore ───────────────────────────────────────────────────────────────
cat > .gitignore << 'EOF'
# Pipeline outputs (large, reproducible)
results/
resources/sequences/
resources/expression/
.snakemake/

# FoldX license (never commit)
foldx_license.txt
rotabase.txt

# Python
__pycache__/
*.pyc
*.egg-info/
.pytest_cache/
dist/

# Conda
*.conda
*.tar.bz2

# OS
.DS_Store
Thumbs.db

# IDE
.vscode/
.idea/
*.swp
EOF
echo "  ✓ .gitignore"

# ── .gitattributes ───────────────────────────────────────────────────────────
cat > .gitattributes << 'EOF'
*.pdb filter=lfs diff=lfs merge=lfs -text
*.cif filter=lfs diff=lfs merge=lfs -text
*.dx  filter=lfs diff=lfs merge=lfs -text
EOF
echo "  ✓ .gitattributes (LFS for structure files)"

# ── LICENSE ──────────────────────────────────────────────────────────────────
cat > LICENSE << 'EOF'
MIT License

Copyright (c) 2026 Pierce Taylor

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOF
echo "  ✓ LICENSE (MIT)"

# ── pyproject.toml ───────────────────────────────────────────────────────────
cat > pyproject.toml << 'EOF'
[project]
name = "carnivorenzyme"
version = "0.1.0"
description = "Structural atlas of convergently evolved digestive enzymes in carnivorous plants"
requires-python = ">=3.11"
license = {text = "MIT"}
authors = [{name = "Pierce Taylor", email = ""}]
dependencies = [
    "biopython>=1.84",
    "pandas>=2.2",
    "numpy>=1.26",
    "scipy>=1.13",
    "matplotlib>=3.9",
    "seaborn>=0.13",
    "click",
    "pyyaml",
    "tqdm",
    "requests",
]

[project.optional-dependencies]
dev = ["pytest", "ruff"]
webapp = ["streamlit>=1.30", "py3Dmol"]
structure = ["chai_lab>=0.6.1"]
docking = ["meeko>=0.5", "vina>=1.2.5"]
evolution = ["evcouplings"]

[tool.ruff]
line-length = 99
target-version = "py311"

[tool.ruff.lint]
select = ["E", "F", "I", "W"]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-v --tb=short"

[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.backends._legacy:_Backend"
EOF
echo "  ✓ pyproject.toml"

# ── environment.yml ──────────────────────────────────────────────────────────
cat > environment.yml << 'YAML'
name: carnivorenzyme
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - snakemake=8.20
  # Bioinformatics
  - mafft=7.525
  - iqtree=2.3.6
  - orthofinder=2.5.5
  - trimal=1.4.1
  - diamond=2.1
  - blast=2.16
  - hmmer=3.4
  # Structure
  - foldseek
  - dssp=4.4.0
  - openbabel=3.1.1
  # Electrostatics
  - apbs=3.4.1
  - pdb2pqr=3.6
  # Expression
  - salmon=1.10.3
  # Python core
  - biopython=1.84
  - pandas=2.2
  - numpy=1.26
  - scipy=1.13
  - matplotlib=3.9
  - seaborn=0.13
  - scikit-learn=1.5
  - click
  - pyyaml
  - tqdm
  - requests
  - pytest
  - ruff
  - pip
  - pip:
    - meeko>=0.5
    - vina>=1.2.5
    - prodigy-prot
    - evcouplings
    - streamlit>=1.30
    - py3Dmol
YAML
echo "  ✓ environment.yml"

# ── config/config.yaml ───────────────────────────────────────────────────────
cat > config/config.yaml << 'YAML'
# CarnivorEnzyme pipeline configuration.
# See CLAUDE.md for full documentation.

project:
  name: CarnivorEnzyme
  version: "0.1.0"

ncbi:
  email: ""     # Set via NCBI_EMAIL env var or directly here
  api_key: ""   # Set via NCBI_API_KEY env var or directly here

orthology:
  inflation: 1.5
  min_species: 3

alignment:
  method: linsi
  trim: automated1

phylogeny:
  model: MFP
  bootstrap: 1000

convergence:
  posterior_threshold: 0.8

structure:
  primary_predictor: alphafold3
  fallback_predictor: chai1
  af3_module: "alphafold3/3.0.1"    # Hellbender module name
  n_seeds: 5
  plddt_exclude: 50

foldx:
  repair_rounds: 5
  temperature: 298
  ddg_destabilizing: 1.0
  ddg_stabilizing: -0.5

evcouplings:
  min_neff: 500           # Minimum effective sequences for reliable DCA
  plmc_iterations: 100

docking:
  exhaustiveness: 32
  num_modes: 20
  box_padding: 10

electrostatics:
  ph_values: [2.5, 3.5, 5.0]
  grid_spacing: 0.5

expression:
  datasets:
    - name: N_gracilis
      sra: PRJDB8591
      genome: GCA_030504385.1
    - name: Cephalotus
      sra: PRJDB4470
      genome: GCA_001941015.1
    - name: Dionaea
      sra: PRJEB12493
      genome: GCA_032404045.1
YAML
echo "  ✓ config/config.yaml"

# ── config/enzyme_families.yaml ──────────────────────────────────────────────
cat > config/enzyme_families.yaml << 'YAML'
# Enzyme families and known accessions.
# Populate accessions before running the pipeline.

tier1:
  chitinases_gh19:
    cazy: GH19
    convergent_source: "Fukushima et al. 2017 Fig. 3a"
    accessions: {}

  purple_acid_phosphatase:
    convergent_source: "Fukushima et al. 2017; Nishimura et al. 2014"
    accessions: {}

  rnase_t2:
    convergent_source: "Fukushima et al. 2017"
    accessions: {}

  nepenthesins:
    merops: A01.073
    subfamily: A1B
    experimental_pdb: null
    accessions:
      Nepenthes_gracilis: ["AB114914", "AB114915"]

  neprosins:
    merops: G03.001
    experimental_pdb: ["7ZVA", "7ZVB", "7ZVC", "7ZU8"]
    accessions: {}

tier2:
  cysteine_proteases: {accessions: {}}
  thaumatin_like: {accessions: {}}
  glucanases_gh17: {accessions: {}}
  lipid_transfer: {accessions: {}}
  esterases: {accessions: {}}
YAML
echo "  ✓ config/enzyme_families.yaml"

# ── SLURM profile ────────────────────────────────────────────────────────────
cat > config/slurm/config.yaml << 'YAML'
# Snakemake SLURM profile for Hellbender HPC.
executor: slurm
default-resources:
  slurm_partition: General
  mem_mb: 8000
  runtime: 120          # minutes
  cpus_per_task: 4
set-resources:
  predict_af3:
    slurm_partition: gpu
    slurm_extra: "'--gres=gpu:1'"
    mem_mb: 64000
    runtime: 480
  run_foldx_scan:
    mem_mb: 4000
    runtime: 60
  run_docking:
    mem_mb: 4000
    runtime: 30
YAML
echo "  ✓ config/slurm/config.yaml"

# ── Snakefile ────────────────────────────────────────────────────────────────
cat > Snakefile << 'PYEOF'
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
PYEOF
echo "  ✓ Snakefile"

# ── Per-rule conda envs ─────────────────────────────────────────────────────
cat > workflow/envs/bioinfo.yaml << 'YAML'
name: bioinfo
channels: [conda-forge, bioconda]
dependencies:
  - mafft=7.525
  - iqtree=2.3.6
  - orthofinder=2.5.5
  - trimal=1.4.1
  - biopython=1.84
  - ete3
  - click
  - pandas
YAML

cat > workflow/envs/structure.yaml << 'YAML'
name: structure
channels: [conda-forge, bioconda]
dependencies:
  - foldseek
  - dssp=4.4.0
  - biopython=1.84
  - click
  - pandas
  - numpy
  - pip
  - pip:
    - chai_lab>=0.6.1
YAML

cat > workflow/envs/foldx.yaml << 'YAML'
name: foldx
channels: [conda-forge]
dependencies:
  - python=3.11
  - pandas
  - click
  - numpy
YAML

cat > workflow/envs/docking.yaml << 'YAML'
name: docking
channels: [conda-forge]
dependencies:
  - openbabel=3.1.1
  - python=3.11
  - click
  - pandas
  - pip
  - pip:
    - meeko>=0.5
    - vina>=1.2.5
YAML

cat > workflow/envs/evcouplings.yaml << 'YAML'
name: evcouplings
channels: [conda-forge, bioconda]
dependencies:
  - python=3.11
  - hmmer=3.4
  - plmc
  - click
  - pandas
  - numpy
  - matplotlib
  - pip
  - pip:
    - evcouplings
YAML

cat > workflow/envs/expression.yaml << 'YAML'
name: expression
channels: [conda-forge, bioconda]
dependencies:
  - salmon=1.10.3
  - click
  - pandas
YAML
echo "  ✓ workflow/envs/"

# ── Placeholder Snakemake rules ──────────────────────────────────────────────
rules=(retrieve orthology alignment phylogeny convergence predict_structure
       structural_align foldx evmutation docking electrostatics expression
       integrate)
for r in "${rules[@]}"; do
  cat > "workflow/rules/${r}.smk" << EOF
# Rule file: ${r}
# TODO: Implement per CLAUDE.md Phase instructions.
EOF
done
echo "  ✓ workflow/rules/ (placeholders)"

# ── Placeholder scripts ──────────────────────────────────────────────────────
scripts=(
  fetch_sequences detect_convergence map_convergence
  predict_chai1 predict_af3 assess_structure classify_positions
  run_foldx_repair run_foldx_scan parse_foldx
  run_evcouplings run_evmutation compare_foldx_evmutation
  prepare_docking run_docking parse_docking
  run_electrostatics quantify_expression
  build_atlas generate_figures
)
for s in "${scripts[@]}"; do
  cat > "workflow/scripts/${s}.py" << PYEOF
#!/usr/bin/env python3
"""${s} — TODO: implement per CLAUDE.md."""

import click


@click.command()
def main():
    raise NotImplementedError("${s} not yet implemented.")


if __name__ == "__main__":
    main()
PYEOF
done
echo "  ✓ workflow/scripts/ (stubs)"

# ── Tests ────────────────────────────────────────────────────────────────────
cat > tests/conftest.py << 'PYEOF'
"""Shared pytest fixtures."""
import pytest
from pathlib import Path

@pytest.fixture
def project_root() -> Path:
    return Path(__file__).parent.parent

@pytest.fixture
def test_data(project_root: Path) -> Path:
    d = project_root / "tests" / "data"
    d.mkdir(exist_ok=True)
    return d
PYEOF

for t in test_fetch test_convergence test_foldx_parse test_docking_parse test_evmutation; do
  cat > "tests/${t}.py" << PYEOF
"""Tests for ${t#test_}. TODO: implement."""

def test_placeholder():
    pass
PYEOF
done
echo "  ✓ tests/"

# ── Webapp stubs ─────────────────────────────────────────────────────────────
cat > webapp/app.py << 'PYEOF'
"""CarnivorEnzyme Streamlit app entry point."""
import streamlit as st

st.set_page_config(page_title="CarnivorEnzyme Atlas", layout="wide")
st.title("CarnivorEnzyme: Structural Atlas")
st.info("Run the Snakemake pipeline first to populate the database.")
PYEOF

for p in family_browser structure_viewer foldx_results docking_results; do
  touch "webapp/pages/${p}.py"
done
echo "  ✓ webapp/"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "=== Repository initialized ==="
echo ""
echo "File count:"
find . -type f | grep -v '.git/' | wc -l
echo ""
echo "Next steps:"
echo "  1. git init && git add -A && git commit -m 'Initial scaffold'"
echo "  2. git remote add origin git@github.com:YOUR_USERNAME/CarnivorEnzyme.git"
echo "  3. git push -u origin main"
echo "  4. conda env create -f environment.yml"
echo "  5. Begin Phase 1: implement workflow/scripts/fetch_sequences.py"
echo ""

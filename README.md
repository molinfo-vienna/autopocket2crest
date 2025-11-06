# AutoPocket2CREST  
> *Automated protein–ligand pocket extraction and preparation pipeline for CREST conformational sampling.*

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Build](https://img.shields.io/badge/status-stable-success.svg)]()
[![CREST](https://img.shields.io/badge/compatible-CREST-orange.svg)](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/crest)

---

## Overview

**AutoPocket2CREST** automates the preparation of **protein–ligand binding pockets** for **CREST** simulations.  
Starting from a protein `.pdb` and ligand `.mol2` file, it automatically:

1. Extracts and cleans the **binding pocket** surrounding the ligand  
2. Adds and optimizes **hydrogen atoms**  
3. Generates **CREST-compatible input files** (`.pdb`, `.xyz`, `.mol2`)  
4. Builds **constraint files** for pocket backbone atoms  
5. Optionally runs **CREST** directly on the prepared system  

The pipeline ensures reproducible, chemically sound local environments ready for quantum mechanical conformer searches.

---

## Features

**Automatic pocket detection** – radius expands dynamically until a sufficient local environment is captured  
**Hydrogenation and protonation** – Open Babel–based hydrogen addition at pH 7.4  
**Structure refinement** – removes alternate locations, unbonded and isolated atoms  
**Constraint generation** – automatically identifies and compresses backbone atom indices for CREST  
**Charge detection** – calculates total formal charge via RDKit  
**CREST-ready output** – `.pdb`, `.xyz`, `.mol2`, and constraints files ready to run  
**Automated cleanup** of intermediate and temporary files  

---

## Input Requirements

| File | Description |
|------|--------------|
| `protein.pdb` | Full protein structure including the bound ligand |
| `ligand.mol2` | Ligand structure file |
| `PDBID` | PDB identifier (used for folder and dataset lookup) |

---

## Generated Output

| File | Purpose |
|------|----------|
| `test_pocket_extended_h_fixed.pdb` | Final CREST-ready pocket structure |
| `test_pocket_extended_h_fixed.xyz` | XYZ file for CREST input |
| `test_pocket_extended_h_fixed.mol2` | MOL2 version for compatibility |
| `constraints.inp` | CREST constraint file |
| `crest.out` | CREST output log |
| `crest_conformers.pdb` | Best conformers in PDB format |

---

## Dependencies

You’ll need the following tools installed and accessible:

| Tool | Purpose |
|------|----------|
| **Python 3.8+** | Core scripting language |
| **MDAnalysis** | Structural analysis |
| **RDKit** | Charge and chemistry operations |
| **Open Babel** | Hydrogenation and file conversion |
| **CREST** | Conformer sampling engine |
| **pdbfixer** | Structure repair and hydrogen addition |
| **OpenMM** | PDB parsing backend |

---

## Installation

```bash
conda create -n autopocket2crest python=3.10
conda activate autopocket2crest
conda install -c conda-forge mdanalysis rdkit openbabel pdbfixer openmm
```

## Usage

Run AutoPocket2CREST directly from the command line:

```bash
python autopocket2crest.py <protein_file.pdb> <ligand_file.mol2> <PDB_ID>
```

Example:

```bash
python autopocket2crest.py 1ABC_protein.pdb 1ABC_ligand.mol2 1ABC
```

This will:

1. Create a working directory /data/local/1ABC

2. Automatically extract and prepare the ligand pocket

3. Generate hydrogenated, CREST-ready structures

4. Run CREST and save results

## Example Workflow

```bash
Input:    1ABC_protein.pdb + 1ABC_ligand.mol2
↓
AutoPocket2CREST
↓
Pocket extraction → Hydrogenation → Cleanup → Constraint generation → CREST
↓
Output: test_pocket_extended_h_fixed.pdb, constraints.inp, crest.out

```

## Why use AutoPocket2CREST

- Fully automated pocket generation and preparation

- Reproducible, standardized quantum-ready systems

- Compatible with CREST workflows

- Simplifies pre-QM preparation for enzyme or receptor-ligand studies

## Author & Citation

Developed by Christian Fellinger (Github: Dragon3221)

This will be part of my Doctoral Thesis. I will share how to Cite it as soon as it is available.

## Acknowledgements

This project builds on tools from the MDAnalysis, RDKit, Open Babel, CREST, and OpenMM communities.
Thank you for all this hard work!
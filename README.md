# ReactionMicroservice
A microservice in Python that would run the reaction SMARTS on reactants supplied in SMILES

## Installation
1. Create a conda environment (skip if you have already an environment with `rdkit`):
```bash
conda env create -f reaction.yml
```
2. Install this package into your environment.
```bash
python3 setup.py install
```

## Usage
This microservice contains one function `run_reaction`, which for given reactants list (SMILES) and reaction SMARTS simulates the reaction:

```python
from ReactionMicroservice import run_reaction

run_reaction(["CC1=CC=C(C=C1)C1=CC(=CC=C1C)C1=CC(C)=CC(C)=C1"], "[c:8]-[c:6]>>[c:8][I:55].[B:99][c:6]")
```
For more details consult docstring:
```python
help(run_reaction)
```

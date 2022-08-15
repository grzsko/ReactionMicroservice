from setuptools import find_packages, setup

setup(
    name='ReactionMicroservice',
    packages=find_packages(),
    version='1.0.0',
    description='A microservice in Python that would run the reaction SMARTS on reactants supplied in SMILES',
    author='Grzegorz Skoraczynski',
    license='MIT',
    python_requires=">=3.9",
    install_requires=["rdkit"]
)

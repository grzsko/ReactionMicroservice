from rdkit import Chem
from rdkit.Chem import AllChem


SMILES = str
SMARTS = str


def run_reaction(reactants: list[SMILES], reaction: SMARTS,
                 products_limit: int = 1) \
        -> list[list[SMILES]]:
    """Run reaction for given list of reactants and reaction.

    Parameters
    ----------
    reactants : iterable of str
        List of SMILES-es representing the reactants. If more than one reactant
        per SMILES, then reactants are separated by comma.
    reaction : str
        SMARTS of reaction tobe run.
    products_limit : positive int
        Maximum number of products sets to be returned
    Returns
    -------
    list of lists of str
        Sets of products. Every set of products contains single product in
        single SMILES.
    """
    # Parse reactants
    mols = []
    try:
        for reactant in reactants:
            reactant = reactant.split(".")
            reactant_mol = [Chem.MolFromSmiles(s) for s in reactant]
            if None in reactant_mol:
                raise ValueError("Something incorrect in SMILES")
            mols.extend([Chem.MolFromSmiles(s) for s in reactant])
    except ValueError as e:
        raise Exception('Incorrect reactant') from e

    # Parse reaction
    try:
        rxn = AllChem.ReactionFromSmarts(reaction)
        if rxn is None:
            raise ValueError("Something incorrect in SMARTS")
    except ValueError as e:
        raise Exception('Incorrect reaction SMART') from e
    
    # Run the reaction
    products = rxn.RunReactants(mols)

    result = []
    for i, product_set in enumerate(products):
        if i == products_limit:
            break
        result.append([Chem.MolToSmiles(s) for s in product_set])

    return result


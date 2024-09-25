"""
A basi fast API for chemistry. 

Implements rdkit and mordred for descriptor generation.
Sam Ellis
"""

from fastapi import FastAPI

import numpt as np

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

from mordred import (
    Calculator,
    descriptors,
    AcidBase,
    HydrogenBond,
    LogS,
    Polarizability,
    SLogP,
    VdwVolumeABC,
)

app = FastAPI()


@app.put("/api/v0/rdkit")
def rdkit_descriptors(smiles: str):
    smiles_clean = sanitise_smiles(smiles)
    mol = smiles_to_mol(smiles_clean)
    cleaned_descriptors = clean_nan_values(descriptors)
    return cleaned_descriptors


@app.put("/api/v0/mordred")
def mordred_descriptors(smiles: str):
    calc = Calculator()
    calc.register(AcidBase.AcidicGroupCount())
    calc.register(AcidBase.BasicGroupCount())
    calc.register(HydrogenBond.HBondAcceptor())
    calc.register(HydrogenBond.HBondDonor())
    calc.register(LogS.LogS())
    calc.register(Polarizability.APol())
    calc.register(Polarizability.BPol())
    calc.register(SLogP.SLogP())
    calc.register(VdwVolumeABC.VdwVolumeABC())
    smiles_clean = sanitise_smiles(smiles)
    mol = smiles_to_mol(smiles_clean)
    descriptors = calc(mol)
    return descriptors.asdict()

# Because some rdkit descriptors are nan, and the fastAPI json_encoder doesn't convert these to null,
# need to do it manually....
def clean_nan_values(data):
    """Recursively replace NaN with None in dict or list."""
    if isinstance(data, dict):
        return {k: clean_nan_values(v) for k, v in data.items()}
    elif isinstance(data, list):
        return [clean_nan_values(item) for item in data]
    elif isinstance(data, float) and np.isnan(data):
        return None  # Replace NaN with None (null in JSON)
    return data


def sanitise_smiles(smiles):
    return smiles


def smiles_to_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol_with_H = Chem.AddHs(mol)
    Chem.AllChem.EmbedMolecule(mol_with_H)
    Chem.AllChem.MMFFOptimizeMolecule(mol_with_H)
    return mol_with_H


# https://greglandrum.github.io/rdkit-blog/posts/2022-12-23-descriptor-tutorial.html
def getMolDescriptors(mol, missingVal=None):
    """calculate the full list of descriptors for a molecule

    missingVal is used if the descriptor cannot be calculated
    """
    res = {}
    for nm, fn in Descriptors._descList:
        # some of the descriptor fucntions can throw errors if they fail, catch those here:
        try:
            val = fn(mol)
        except:
            # print the error message:
            import traceback

            traceback.print_exc()
            # and set the descriptor value to whatever missingVal is
            val = missingVal
        res[nm] = val
    return res

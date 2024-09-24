"""
A basic flask API for chemistry. 

Implements rdkit and mordred for descriptor generation.
Sam Ellis
"""

from flask import Flask, request, jsonify

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

import numpy as np

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

app = Flask(__name__)


@app.route("/api/v0/descriptors/rdkit", methods=["POST"])
def rdkitdescriptors():
    data = request.get_json()
    if data["smiles"] is not None:
        try:
            smiles = sanitise_smiles(data["smiles"])
            mol = smiles_to_mol(smiles)

            results = getMolDescriptors(mol)

            return (
                jsonify(results),
                200,
            )
        except Exception as e:
            return (
                jsonify(
                    {
                        "error": str(e),
                    }
                ),
                400,
            )

    else:
        return jsonify({"error": "Invalid JSON data"}), 400


@app.route("/api/v0/mordred", methods=["POST"])
def mordred_descriptors():
    data = request.get_json()
    if data["smiles"] is not None:
        try:
            smiles = sanitise_smiles(data["smiles"])
            mol = smiles_to_mol(smiles)
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
            mordred_results = calc(mol)
            return (
                jsonify(mordred_results.asdict()),
                200,
            )
        except Exception as e:
            return (
                jsonify(
                    {
                        "error": str(e),
                    }
                ),
                400,
            )

    else:
        return jsonify({"error": "Invalid JSON data"}), 400


# Get a specific mordred descriptor
@app.route("/api/v0/descriptors/mordred/<desc>", methods=["POST"])
def mordred_descriptors_single(desc):
    mordred_descriptors = {
        "AcidicGroupCount": AcidBase.AcidicGroupCount(),
        "BasicGroupCount": AcidBase.BasicGroupCount(),
    }
    data = request.get_json()
    if data["smiles"] is not None:
        try:
            smiles = sanitise_smiles(data["smiles"])
            mol = smiles_to_mol(smiles)
            calc = Calculator()
            calc.register(mordred_descriptors.get(desc))
            mordred_results = calc(mol)

            results = mordred_results.asdict()
            return (
                jsonify(results),
                200,
            )
        except Exception as e:
            return (
                jsonify(
                    {
                        "error": str(e),
                    }
                ),
                400,
            )

    else:
        return jsonify({"error": "Invalid JSON data"}), 400


@app.route("/api/v0/descriptors/morgan", methods=["POST"])
def morgan():
    data = request.get_json()
    if data["smiles"] is not None:
        try:
            smiles = sanitise_smiles(data["smiles"])
            mol = smiles_to_mol(smiles)
            mfpt = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=12)
            mfpt_arr = np.array(mfpt)
            return (
                jsonify(
                    {
                        "results": str(mfpt_arr),
                    }
                ),
                200,
            )
        except Exception as e:
            return (
                jsonify(
                    {
                        "error": str(e),
                    }
                ),
                400,
            )


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


def sanitise_smiles(smiles):
    return smiles


def smiles_to_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol_with_H = Chem.AddHs(mol)
    Chem.AllChem.EmbedMolecule(mol_with_H)
    Chem.AllChem.MMFFOptimizeMolecule(mol_with_H)
    return mol_with_H


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)

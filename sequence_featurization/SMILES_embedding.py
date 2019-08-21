import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import pairwise_distances
import warnings

drugs_filename = "../data/compounds.txt"
reference_drugs_filename = "reference_SMILES.csv"

size_limit = 15000


def get_fingerprints(smiles_string):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    except:
        return np.nan


def main():
    dataset_smiles = pd.read_csv(drugs_filename).set_index("cid").smiles
    dataset_fingerprints = np.array(list(dataset_smiles.apply(get_fingerprints)))
    reference_smiles = pd.read_csv(reference_drugs_filename, header=-1).set_index(0)[1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reference_fingerprints = reference_smiles.apply(get_fingerprints)
        reference_fingerprints = np.array(list(reference_fingerprints.dropna()))
        reference_fingerprints = pd.DataFrame(reference_fingerprints).drop_duplicates().sample(size_limit, random_state=42)
        reference_distances = pairwise_distances(np.vstack([reference_fingerprints, dataset_fingerprints]),
                                                 metric="jaccard", n_jobs=-1)[:, :size_limit]
    pd.DataFrame(reference_distances).to_csv("reference_SMILES_distances.csv")


if __name__ == '__main__':
    main()

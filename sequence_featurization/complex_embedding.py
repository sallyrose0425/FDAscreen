import numpy as np
import pandas as pd
import multiprocessing as mp
from skbio.alignment import StripedSmithWaterman
from umap import UMAP
from Bio.SubsMat import MatrixInfo
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import pairwise_distances


def subs_mat(subs):
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                   'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X']
    a_a_pairs = {}
    for c in amino_acids:
        a_a_pairs[c] = {'*': 0}
        for s in amino_acids:
            try:
                a_a_pairs[c][s] = subs[(c, s)]
            except KeyError:
                a_a_pairs[c][s] = subs[(s, c)]
    dummy_dict = {'*': {'*': 0}}
    for c in amino_acids:
        dummy_dict['*'][c] = 0
    a_a_pairs.update(dummy_dict)
    return a_a_pairs


def fasta_metric(s, S):
    query = StripedSmithWaterman(s, protein=True, substitution_matrix=subs_mat(MatrixInfo.pam250))
    scores = []
    for t in S:
        scores.append(query(t)['optimal_alignment_score'])
    return scores

"""
def get_fingerprints(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
"""

complexes_filename = "sequence_featurization/complexes.csv"
complexes = pd.read_csv(complexes_filename)
prot_unique_mask = complexes.prot_ID.unique()
complexes[prot_unique_mask]
prots_ids = prots.tolist()[:20]


def distances_list(index):
    distance_list = fasta_metric(prots[index], prots)
    return [dist / distance_list[index] for dist in distance_list]


def main():
    pool = mp.Pool()
    distances = pool.map(distances_list, range(len(prots)))
    fasta_distances = 1 - np.array(distances)
    umap = UMAP(n_components=10)
    fasta_embedding = umap.fit_transform(fasta_distances)



""" # load data
fastas = pd.read_csv('data/proteins.txt', index_col=0).set_index('pid')
smiles = pd.read_csv('data/compounds.txt', index_col=0).set_index('cid')
interactions = pd.read_csv('data/interactions.txt', index_col=0)


# generate FASTA distance matrix and embedding
fasta_list = fastas.fasta.tolist()



# generate fingerprint distance matrix and embedding
fingerprints = pd.DataFrame(np.array(smiles.smiles.apply(get_fingerprints).tolist()))
fingerprint_distances = pairwise_distances(fingerprints.values.astype(bool), metric='jaccard')
"""

if __name__ == '__main__':
    main()
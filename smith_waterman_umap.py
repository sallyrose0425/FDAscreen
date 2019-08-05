import numpy as np
import pandas as pd
from skbio.alignment import StripedSmithWaterman
from umap import UMAP
from Bio.SubsMat import MatrixInfo
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import Normalizer


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


def get_fingerprints(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


# load data
fastas = pd.read_csv('data/proteins.txt', index_col=0).set_index('pid')
smiles = pd.read_csv('data/compounds.txt', index_col=0).set_index('cid')
interactions = pd.read_csv('data/interactions.txt', index_col=0)


# generate FASTA distance matrix and embedding
fasta_list = fastas.fasta.tolist()
distance_lists = []
for index in range(len(fasta_list)):
    distance_list = fasta_metric(fasta_list[index], fasta_list)
    distance_lists.append([dist / distance_list[index] for dist in distance_list])
fasta_distances = 1 - np.array(distance_lists)
fasta_embedding = Normalizer().fit_transform(UMAP(n_components=10).fit_transform(fasta_distances))
fasta_embedding = pd.DataFrame(fasta_embedding).set_index(fastas.index)


# generate fingerprint distance matrix and embedding
fingerprints = pd.DataFrame(np.array(smiles.smiles.apply(get_fingerprints).tolist()))
fingerprint_distances = pairwise_distances(fingerprints.values.astype(bool), metric='jaccard')
fingerprint_embedding = Normalizer().fit_transform(UMAP(n_components=10).fit_transform(fingerprint_distances))
fingerprint_embedding = pd.DataFrame(fingerprint_embedding).set_index(smiles.index)


# generate data feature from combined distances
def combined_features(row):
    return fingerprint_embedding.loc[row.cid].tolist() + fasta_embedding.loc[row.pid].tolist()


data_features = pd.DataFrame(interactions.apply(combined_features, axis=1).tolist())


# visualization
activity_filter = interactions.activity.astype(bool).values
data_viz = UMAP().fit_transform(data_features)
plt.scatter(data_viz[:, 0][~activity_filter], data_viz[:, 1][~activity_filter], alpha=0.05)
plt.scatter(data_viz[:, 0][activity_filter], data_viz[:, 1][activity_filter], alpha=0.05)


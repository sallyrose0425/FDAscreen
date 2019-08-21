import numpy as np
import pandas as pd
import multiprocessing as mp
from skbio.alignment import StripedSmithWaterman
from umap import UMAP
from Bio.SubsMat import MatrixInfo
import time

size_limit = 5000
reference_fastas = pd.read_csv("reference_FASTAS.csv", header=-1,
                               names=["ID", "FASTA"]).set_index("ID").sample(n=size_limit, random_state=42)
dataset_fastas = pd.read_csv("../data/proteins.txt", header=0, names=["ID", "FASTA"]).set_index("ID")
fastas = pd.concat([reference_fastas, dataset_fastas], sort=False)


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


substitution_matrix = subs_mat(MatrixInfo.pam250)


def fasta_reference_similarity(s):
    query = StripedSmithWaterman(s, protein=True, substitution_matrix=substitution_matrix)
    return reference_fastas.apply(lambda t: query(t.FASTA)['optimal_alignment_score'], axis=1).values


def main():
    t0 = time.time()
    pool = mp.Pool()
    similarities = np.array(pool.map(fasta_reference_similarity, fastas.FASTA.values.tolist()))
    norms = np.diagonal(similarities[:size_limit, :])
    reference_distances = 1 - (similarities / norms[:, None].transpose())
    reference_distances = pd.DataFrame(reference_distances).set_index(fastas.index)
    reference_distances.to_csv("reference_distances.csv")
    print(time.time() - t0)


if __name__ == '__main__':
    main()

import pandas as pd

fastas = pd.read_csv("sequence_featurization/FASTAS.csv", header=None, names=["prot_ID", "FASTA"])
smiles = pd.read_csv("sequence_featurization/SMILES.csv", header=None, names=["lig_ID", "SMILES"])

size_cap = 15000  # too many samples and the distance matrix won't fit in memory...
smiles_samples = smiles.sample(n=size_cap, random_state=42).reset_index().drop("index", axis=1)
fastas_samples = pd.concat([fastas]*((size_cap // len(fastas)) + 1)).reset_index().drop("index", axis=1)[:size_cap]
complexes = smiles_samples.join(fastas_samples).set_index(["prot_ID", "lig_ID"])
complexes.to_csv("sequence_featurization/complexes.csv")

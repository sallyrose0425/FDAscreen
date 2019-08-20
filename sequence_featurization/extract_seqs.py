import pandas as pd

#####################################################################
fastas = {}
with open("sequence_featurization/UP000005640_9606.fasta", "r") as fastas_file:
    text = fastas_file.read()
prots = iter(text.split(">")[1:])
for prot in prots:
    try:
        ID = prot.split(" ")[0]
        FASTA = prot[prot.index("\n"):].replace("\n", "")
        fastas[ID] = FASTA
    except ValueError:
        continue

fastas_frame = pd.Series(fastas)
fastas_frame.to_csv("sequence_featurization/reference_FASTAS.csv", header=False)

#####################################################################
smiles = {}
with open("sequence_featurization/structures.sdf", "r") as smiles_file:
    text = smiles_file.read()
mols = iter(text.split("> <DATABASE_ID>")[1:])
for mol in mols:
    ligand = mol.split("\n")[1]
    SMILES = mol.split("<SMILES>\n")[1].split("\n")[0]
    smiles[ligand] = SMILES
smiles_frame = pd.Series(smiles)
smiles_frame.to_csv("sequence_featurization/reference_SMILES.csv", header=False)

import pandas as pd


#####################################################################
fastas = {}
target = None
with open("sequence_featurization/protein.fasta", "r") as fastas_file:
    for line in fastas_file.readlines():
        if line[0] == ">":
            target = line.split(" ")[0][1:]
            fastas[target] = ""
        else:
            line = line.replace("\n", "")
            fastas[target] += line
fastas_frame = pd.Series(fastas)
fastas_frame.to_csv("sequence_featurization/FASTAS.csv", header=False)

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
smiles_frame.to_csv("sequence_featurization/SMILES.csv", header=False)


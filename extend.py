from __future__ import print_function
import csv
from util import *

datadir = 'data/consolidated/'
#datadir = 'data/'
drugs,proteins,interactions = {}, {}, {}
dr = csv.reader(open(datadir + 'compounds.txt', 'r'))
pr = csv.reader(open('data/extended/' + 'proteins.txt', 'r'))
ir = csv.reader(open(datadir + 'interactions.txt', 'r'))
print(next(dr), next(pr), next(ir))
for i, cid, smiles in dr:
    drugs[cid] = smiles
for i, pid, sequence in pr:
    proteins[pid] = sequence
for i, cid, pid, activity in ir:
    interactions[(cid, pid)] = int(activity)
drugs_,proteins_,interactions_ = map(dict, (drugs, proteins, interactions))

'''
new_interactions = {}
bad_pids = []
# Query proteins for all known drugs, adding proteins as they are found.
for cid in drugs:
    table = pubchem_query_cid_interactions(cid)
    if 'protacxn' in table:
        for row in table[['protacxn', 'activity']].itertuples():
            if row.protacxn == 'NULL' or row.protacxn in bad_pids or row.activity not in ('Active', 'Inactive'):
                continue
            protacxn = row.protacxn
            if protacxn not in proteins:
                print(cid,protacxn)
                try:
                    sequence = pubchem_pid_to_fasta(protacxn)
                except Exception as e:
                    print(e)
                    bad_pids.append(protacxn)
                    continue
                proteins[protacxn] = sequence
            else:
                sequence = proteins[row.protacxn]
            if (cid,pid) in new_interactions:
                actives, inactives = new_interactions[(cid,pid)]
            else:
                actives, inactives = 0, 0
            if row.activity == 'Active':
                actives += 1
            elif row.activity == 'Inactive':
                inactives += 1
            new_interactions[(cid,pid)] = actives, inactives

for ((cid, pid), (actives, inactives)) in new_interactions.items():
    if (cid, pid) in interactions and interactions[(cid,pid)] == 1:
        activity = 1
    elif actives < inactives:
        activity = 0
    else:
        activity = 1
    interactions[(cid,pid)] = activity
'''

datadir = 'data/extended/'

with open(datadir + 'proteins.txt', 'w') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(('', 'pid', 'sequence'))
    for index, (pid, sequence) in enumerate(proteins.items()):
        w.writerow((index, pid, sequence))

with open(datadir + 'interactions.txt', 'w') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(('', 'cid', 'pid', 'activity'))
    for index, ((cid, pid), activity) in enumerate(interactions.items()):
        w.writerow((index, cid, pid, activity))

new_interactions = {}
pids = list(proteins.keys())
# Use already queried drug info to save time.
dr = csv.reader(open(datadir + 'compounds.txt', 'r'))
next(dr)
for i, cid, smiles in dr:
    drugs[cid] = smiles
# Query interactions for all proteins, adding drugs as they are found.
new_count = 0
while pids:
    pid = pids.pop(0)
    try:
        rows = pubchem_query_protein_interactions(pid)
    except AssertionError as e:
        print('Error querying protein interactions', str(e))
        pids.append(pid)
        continue
    except json.JSONDecodeError as e:
        print('JSONDecodeError', str(e))
        pids.append(pid)
        continue
    for row in rows:
        cid, activity = row['cid'], row['activity']
        if activity not in ('Active', 'Inactive'):
            continue
        activity = int(activity == 'Active')
        if cid not in drugs:
            try:
                smiles = pubchem_cid_to_canonical_smiles(cid)
                drugs[cid] = smiles
                new_count += 1
            except AssertionError as e:
                print('Error getting SMILES for cid', cid)
                continue
        else:
            smiles = drugs[cid]
        if (cid, pid) in new_interactions:
            inactives, actives = new_interactions[(cid,pid)]
        else:
            inactives, actives = 0, 0
        if activity == 0: inactives += 1
        if activity == 1: actives += 1
        new_interactions[(cid,pid)] = (inactives, actives)
        sequence = proteins[pid]
        print(pid, cid, actives, inactives)
        if new_count > 1000:
            with open(datadir + 'compounds.txt', 'w') as f:
                w = csv.writer(f, lineterminator='\n')
                w.writerow(('', 'cid', 'smiles'))
                for index, (cid, smiles) in enumerate(drugs.items()):
                    w.writerow((index, cid, smiles))
            new_count = 0


for ((cid, pid), (actives, inactives)) in new_interactions.items():
    if (cid, pid) in interactions and interactions[(cid,pid)] == 1:
        activity = 1
    elif actives < inactives:
        activity = 0
    else:
        activity = 1
    interactions[(cid,pid)] = activity

with open(datadir + 'compounds.txt', 'w') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(('', 'cid', 'smiles'))
    for index, (cid, smiles) in enumerate(drugs.items()):
        w.writerow((index, cid, smiles))

with open(datadir + 'proteins.txt', 'w') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(('', 'pid', 'sequence'))
    for index, (pid, sequence) in enumerate(proteins.items()):
        w.writerow((index, pid, sequence))

with open(datadir + 'interactions.txt', 'w') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(('', 'cid', 'pid', 'activity'))
    for index, ((cid, pid), activity) in enumerate(interactions.items()):
        w.writerow((index, cid, pid, activity))

from __future__ import print_function
import csv

#datadir = 'data/consolidated/'
datadir = 'data/extended/'
interactions = {}
with open(datadir + 'interactions.txt', 'r') as f:
    r = csv.reader(f)
    next(r)
    for i, cid, pid, activity in r:
        interactions[(cid, pid)] = int(activity)

drugs, proteins = {}, {}
active, total = 0, 0
for (cid, pid), activity in interactions.items():
    if cid not in drugs:
        drugs[cid] = (0, 0)
    if pid not in proteins:
        proteins[pid] = (0, 0)
    drug_active, drug_total = drugs[cid]
    protein_active, protein_total = proteins[pid]
    if activity:
        drug_active += 1
        protein_active += 1
        active += 1
    drug_total += 1
    protein_total += 1
    drugs[cid] = (drug_active, drug_total)
    proteins[pid] = (protein_active, protein_total)
    total += 1
print('Active: %d(%d%%)\n\nTotal: %d\n\n' % (active, float(active)/total*100, total))

print('\nDrugs\n')
only_active, only_inactive, few_interactions = 0, 0, 0
for cid in drugs:
    active, total = drugs[cid]
    if active == total:
        only_active += 1
        #print('Only Active: ', cid)
    if active == 0:
        only_inactive += 1
    if total < 5:
        few_interactions += 1
    #print('%s: %d/%d (%d%%)\n' % (cid, active, total, float(active)/total*100))
total = len(drugs)
print('Total: %d\n\nOnly Actives: %d(%d%%)\n\nOnly Inactives: %d(%d%%)\n\nFewer than 5 interactions: %d(%d%%)\n\n' % 
      (total,
       only_active, float(only_active)/total*100, 
       only_inactive, float(only_inactive)/total*100, 
       few_interactions, float(few_interactions)/total*100))

print('\nProteins\n')
only_active, only_inactive, few_interactions = 0, 0, 0
for pid in proteins:
    active, total = proteins[pid]
    if active == total:
        only_active += 1
    if active == 0:
        only_inactive += 1
    if total < 5:
        few_interactions += 1
    #print('%s: %d/%d (%d%%)\n' % (pid, active, total, float(active)/total*100))
total = len(proteins)
print('Total: %d\n\nOnly Actives: %d(%d%%)\n\nOnly Inactives: %d(%d%%)\n\nFewer than 5 interactions: %d(%d%%)\n\n' % 
      (total,
       only_active, float(only_active)/total*100, 
       only_inactive, float(only_inactive)/total*100, 
       few_interactions, float(few_interactions)/total*100))

import csv

interactions = {}
with open('data/consolidated/interactions.txt', 'r') as f:
    r = csv.reader(f)
    r.next()
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
print('Active: %d(%d%%)\n\n Total: %d\n\n' % (active, float(active)/total*100, total))
print('\nDrugs\n')
for cid in drugs:
    active, total = drugs[cid]
    print('%s: %d/%d (%d%%)\n' % (cid, active, total, float(active)/total*100))
print('\nProteins\n')
for pid in proteins:
    active, total = proteins[pid]
    print('%s: %d/%d (%d%%)\n' % (pid, active, total, float(active)/total*100))

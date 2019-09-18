from __future__ import print_function
import os, csv, requests
from util import pubchem_smiles_to_cid, read_broken_file, strip_fasta

# https://www.ebi.ac.uk/Tools/sss/fasta/
# A possible search tool for FASTA.
# https://blast.ncbi.nlm.nih.gov/Blast.cgi
# This is giving completely different results.

def extract_ttd():
    # A dictionary mapping TTD Drug ID to (Drug Name, PubChem CID, SMILES, InChi).
    ttd_drugs = {}
    # A dictionary mapping TTD Target ID to (UniProt ID, Target Name, Sequence).
    ttd_targets = {}
    # A dictionary mapping (TTD Drug ID, TTD Target ID) (Target Name, Actions)
    ttd_interactions = {}
    # A dictionary mapping Pubchem CID to TTD Drug ID.
    drug_pubchem_to_ttd = {}
    # A dictionary mapping UniProt PID to TTD Target ID
    protein_uniprot_to_ttd = {}

    ### bidd.nus ###
    ### bidd.nus/drugDataset.txt ###
    # Files do have column headers (need to skip first line).
    # Three of the entries in drugDataset.txt don't have PubChemCIDs.
    #  pubchem_smiles_to_cid is able to find PubChemCIDs for all three of them.
    print('Processing bidd.nus/drugDataset.txt...')
    drugs = csv.reader(open('data/bidd.nus/drugDataset.txt', 'r'), delimiter='\t')
    next(drugs) # Skip column header line.
    lines, good_lines = 0, 0 # Keep track of how many entries successfully processed.
    for ttd_drugid, drugname, pubchem_cid, smiles, inchi in drugs:
        lines += 1
        assert ttd_drugid not in ttd_drugs # Assert that there are no duplicate entries.
        assert smiles # Assert that there is a SMILES string.
        if not pubchem_cid: # Look up PubChemCID based on SMILES if it's not there.
            warning = ('Missing PubChem CID', (ttd_drugid, drugname, pubchem_cid, smiles, inchi))
            pubchem_cid = pubchem_smiles_to_cid(smiles)
        assert pubchem_cid
        try:
            assert pubchem_cid not in drug_pubchem_to_ttd
        except AssertionError as e:
            old_drugid = drug_pubchem_to_ttd[pubchem_cid]
            old = (old_drugid,) + ttd_drugs[old_drugid]
            old_smiles = old[3]
            if smiles != old_smiles:
                error = ('Duplicate PubChem CID with non-matching SMILES', 
                         old, (ttd_drugid, drugname, pubchem_cid, smiles, inchi))
                #print(repr(error))
                continue
            else:
                # The new SMILES matches the old SMILES, include in TTD Drugs data and print a warning.
                warning = ('Duplicate PubChem CID with matching SMILES', 
                           old, (ttd_drugid, drugname, pubchem_cid, smiles, inchi))
                #print(repr(warning))
        ttd_drugs[ttd_drugid] = (drugname, pubchem_cid, smiles, inchi)
        drug_pubchem_to_ttd[pubchem_cid] = ttd_drugid
        good_lines += 1
    print('Processed %d of %d lines' % (good_lines, lines))
    ### bidd.nus/targetDataset.txt ###
    # Many of the entries have residue numbers.
    # One entry does not have a UniProtID at all.
    #  Manual searches have not found a matching FASTA.
    # The sequences are listed without FASTA header.
    print('Processing bidd.nus/targetDataset.txt...')
    targets = csv.reader(open('data/bidd.nus/targetDataset.txt', 'r'), delimiter='\t')
    next(targets) # Skip column header line.
    lines, good_lines = 0, 0 # Keep track of how many entries successfully processed.
    for ttd_targetid, uniprotid, targetname, sequence in targets:
        lines += 1
        assert ttd_targetid not in ttd_targets
        try:
            assert uniprotid
        except AssertionError as e:
            error = ('Missing UniProt ID', (ttd_targetid, uniprotid, targetname, sequence))
            print(repr(error))
            continue
        assert sequence
        assert ttd_targetid not in ttd_targets
        try:
            assert uniprotid not in protein_uniprot_to_ttd
        except AssertionError as e:
            old_targetid = protein_uniprot_to_ttd[uniprotid]
            old = (old_targetid,) + ttd_targets[old_targetid]
            old_sequence = old[3]
            if sequence != old_sequence:
                error = ('Duplicate UniProt ID with non-matching sequence', 
                         old, (ttd_targetid, uniprotid, targetname, sequence))
                #print(repr(error))
                continue
            else:
                warning = ('Duplicate UniProt ID with matching sequence', 
                           old, (ttd_targetid, uniprotid, targetname, sequence))
                #print(repr(warning))
        ttd_targets[ttd_targetid] = (uniprotid, targetname, sequence)
        protein_uniprot_to_ttd[uniprotid] = ttd_targetid
        good_lines += 1
    print('Processed %d of %d lines' % (good_lines, lines))
    ### bidd.nus/interactionDataset.txt ###
    print('Processing bidd.nus/interactionDataset.txt...')
    interactions = csv.reader(open('data/bidd.nus/interactionDataset.txt', 'r'), delimiter='\t')
    next(interactions) # Skip column header line.
    lines, good_lines = 0, 0 # Keep track of how many entries successfully processed.
    for ttd_drugid, ttd_targetid, targetname, actions in interactions:
        lines += 1
        try:
            assert ttd_drugid in ttd_drugs
        except AssertionError as e:
            error = ('Unknown Drug ID', ttd_drugid)
            print(repr(error))
            continue
        try:
            assert ttd_targetid in ttd_targets
        except AssertionError as e:
            error = ('Unknown Target ID', ttd_targetid)
            print(repr(error))
            continue
        assert (ttd_drugid, ttd_targetid) not in ttd_interactions
        ttd_interactions[(ttd_drugid, ttd_targetid)] = (targetname, actions)
        good_lines += 1
    print('Processed %d of %d lines' % (good_lines, lines))
    return ttd_drugs, ttd_targets, ttd_interactions

# Takes the dictionaries extracted by extract_ttd and
# puts into format suitable for combining with other data.
def process_ttd(ttd_drugs, ttd_targets, ttd_interactions):
    drugs, proteins, interactions = [], [], []
    for ttd_drugid, (drug_name, pubchem_cid, smiles, inchi) in ttd_drugs.items():
        assert pubchem_cid
        assert smiles
        drugs.append((pubchem_cid, smiles))
    for ttdtargetid, (uniprot_id, targetname, sequence) in ttd_targets.items():
        assert uniprot_id
        assert sequence
        proteins.append((uniprot_id, sequence))
    for (ttd_drugid, ttd_targetid), (targetname, actions) in ttd_interactions.items():
        drug_name, pubchem_cid, smiles, inchi = ttd_drugs[ttd_drugid]
        uniprot_id, target_name, sequence = ttd_targets[ttd_targetid]
        interactions.append((pubchem_cid, uniprot_id, 1))
    return drugs, proteins, interactions

def extract_db():
    # A dictionary mapping DB Drug ID to (Drug Name, SMILES, InChi, PubChem CID).
    db_drugs = {}
    # A dictionary mapping DB Protein ID to (UniProt ID, UniProtKB Entry Name, Protein Name, FASTA, Sequence).
    db_proteins = {}
    # A dictionary mapping Sequence to (DB Protein ID, Protein Name).
    db_proteins_by_sequence = {}
    # A dictionary mapping (TTD Drug ID, TTD Target ID) to (Target Name, Actions)
    db_interactions = {}

    ### drugbankversionmonth12Final ###
    # Files do not have column headers.
    # Lines in proteinDataset.txt contain newlines, breaking usual tab-delimited file reading code.
    # I've built a special routine for reading that file.
    print('Processing drugbankversionmonth12Final/drugDataset.txt...')
    drugs = csv.reader(open('data/drugbankversionmonth12Final/drugDataset.txt', 'r'), delimiter='\t')
    lines, good_lines = 0, 0 # Keep track of how many entries successfully processed.
    for db_drugid, drug_name, smiles, inchi, pubchem_cid in drugs:
        lines += 1
        assert db_drugid not in db_drugs
        assert smiles
        if not pubchem_cid: # Look up PubChemCID based on SMILES if it's not there.
            warning = ('Missing PubChem CID', (db_drugid, drug_name, smiles, inchi, pubchem_cid))
            print(repr(warning))
            try:
                # Fails on 'Cl[223Ra]Cl'
                pubchem_cid = pubchem_smiles_to_cid(smiles)
            except Exception as e:
                print(e)
                continue
            continue
        db_drugs[db_drugid] = (drug_name, smiles, inchi, pubchem_cid)
        good_lines += 1
    print('Processed %d of %d lines' % (good_lines, lines))
    1/0

    print('Processing drugbankversionmonth12Final/proteinDataset.txt...')
    proteins = read_broken_file('data/drugbankversionmonth12Final/proteinDataset.txt')
    count = 0
    for db_proteinid, uniprot_id, uniprotkb_entryname, protein_name, fasta in proteins:
        try:
            assert (db_proteinid, protein_name) not in db_proteins
        except AssertionError as e:
            print('Duplicate protein entry:')
            print((db_proteinid, protein_name), (uniprot_id, uniprotkb_entryname, fasta))
            print((db_proteinid, protein_name), db_proteins[(db_proteinid, protein_name)][:-1])
            continue
        assert uniprot_id
        assert fasta
        sequence = strip_fasta(fasta)
        db_proteins[(db_proteinid, protein_name)] = (uniprot_id, uniprotkb_entryname, fasta, sequence)
        try:
            assert sequence not in db_proteins_by_sequence
        except AssertionError as e:
            print('Duplicate sequence', sequence)
            print(db_proteins_by_sequence[sequence])
            print((db_proteinid, protein_name))
        db_proteins_by_sequence[sequence] = (db_proteinid, protein_name)
        count += 1
    print('Processed %d lines' % count)

    print('Processing drugbankversionmonth12Final/interactionsDataset.txt...')
    interactions = csv.reader(open('data/drugbankversionmonth12Final/interactionsDataset.txt', 'r'), delimiter='\t')
    count = 0
    for db_drugid, db_proteinid, receptor, pharma_action, protein_name, actions in interactions:
        try:
            assert db_drugid in db_drugs
        except AssertionError as e:
            print('Missing db_drugid', db_drugid)
            continue
        try:
            assert (db_proteinid, protein_name) in db_proteins
        except AssertionError as e:
            print('Missing (db_proteinid, protein_name)', db_proteinid, protein_name)
            continue
        try:
            assert (db_drugid, (db_proteinid, protein_name)) not in db_interactions
        except AssertionError as e:
            print('Duplicate interaction')
            print((db_drugid, (db_proteinid, protein_name)), (receptor, pharma_action, actions))
            print((db_drugid, (db_proteinid, protein_name)), db_interactions[(db_drugid, (db_proteinid, protein_name))])
            print('*'*20)
        db_interactions[(db_drugid, (db_proteinid, protein_name))] = (receptor, pharma_action, actions)
        count += 1
    print('Processed %d lines' % count)

    return db_drugs, db_proteins, db_interactions

# Takes the dictionaries extracted by extract_ttd and
# puts into format suitable for combining with other data.
def process_db(db_drugs, db_proteins, db_interactions):
    drugs, proteins, interactions = [], [], []
    for db_drugid, (drug_name, smiles, inchi, pubchem_cid) in db_drugs.items():
        assert pubchem_cid
        assert smiles
        drugs.append((pubchem_cid, smiles))
    for (db_proteinid, protein_name), (uniprot_id, uniprotkb_entryname, fasta, sequence) in db_proteins.items():
        assert uniprot_id
        assert sequence
        proteins.append((uniprot_id, sequence))
    for (db_drugid, (db_proteinid, protein_name)), (receptor, pharma_action, actions) in db_interactions.items():
        drug_name, smiles, inchi, pubchem_cid = db_drugs[db_drugid]
        uniprot_id, uniprotkb_entryname, fasta, sequence = db_proteins[(db_proteinid, protein_name)]
        interactions.append((pubchem_cid, uniprot_id, 1))
    return drugs, proteins, interactions


#TODO: Add checking for duplicate smiles.
def merge_drugs(rows1, rows2):
    drugs = {}
    for cid, smiles in rows1:
        assert cid not in drugs
        drugs[cid] = smiles
    for cid, smiles in rows2:
        if cid in drugs:
            try:
                assert drugs[cid] == smiles
            except AssertionError as e:
                print('DUPLICATE CID WITH NONMATCHING SMILES', cid, drugs[cid], smiles)
                continue
            continue
        drugs[cid] = smiles
    return drugs

def merge_proteins(rows1, rows2):
    proteins = {}
    for pid, sequence in rows1:
        assert pid not in proteins
        proteins[pid] = sequence
    for pid, sequence in rows2:
        if pid in proteins:
            try:
                assert proteins[pid] == sequence
            except AssertionError as e:
                print('DUPLICATE PID WITH NONMATCHING SEQUENCE', pid, sequence, proteins[pid])
            continue
        proteins[pid] = sequence
    return proteins

def merge_interactions(rows1, rows2):
    interactions = {}
    for cid, pid, activity in rows1:
        assert (cid, pid) not in interactions
        interactions[(cid, pid)] = activity
    for cid, pid, activity in rows2:
        if (cid, pid) in interactions:
            try:
                assert interactions[(cid, pid)] == activity
            except AssertionError as e:
                print('Conflicting activity data for', (cid, pid))
            continue
        interactions[(cid, pid)] = activity
    return interactions

def main():
    drugs = []
    with open('data/compounds.txt', 'r') as f:
        r = csv.reader(f)
        next(r)
        for index, cid, smiles in r:
            drugs.append((cid, smiles))

    proteins = []
    with open('data/proteins.txt', 'r') as f:
        r = csv.reader(f)
        next(r)
        for index, pid, sequence in r:
            proteins.append((pid, sequence))

    interactions = []
    with open('data/interactions.txt', 'r') as f:
        r = csv.reader(f)
        next(r)
        for index, cid, pid, activity in r:
            interactions.append((cid, pid, activity))

    ttd_drugs, ttd_targets, ttd_interactions = extract_ttd()
    db_drugs, db_proteins, db_interactions = extract_db()
    drugs1, proteins1, interactions1 = process_ttd(ttd_drugs, ttd_targets, ttd_interactions)
    drugs2, proteins2, interactions2 = process_db(db_drugs, db_proteins, db_interactions)
    drugs = merge_drugs(drugs, drugs1)
    drugs = merge_drugs(drugs.items(), drugs2)
    proteins = merge_proteins(proteins, proteins1)
    proteins = merge_proteins(proteins.items(), proteins2)
    interactions = merge_interactions(interactions, interactions1)
    interactions = merge_interactions([(cid, pid, activity) for (cid, pid), activity in interactions.items()], interactions2)

    with open('data/consolidated/compounds.txt', 'w') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(('', 'cid', 'smiles'))
        for index, (cid, smiles) in enumerate(drugs.items()):
            w.writerow((index, cid, smiles))

    with open('data/consolidated/proteins.txt', 'w') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(('', 'pid', 'sequence'))
        for index, (pid, sequence) in enumerate(proteins.items()):
            w.writerow((index, pid, sequence))

    with open('data/consolidated/interactions.txt', 'w') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(('', 'cid', 'pid', 'activity'))
        for index, ((cid, pid), activity) in enumerate(interactions.items()):
            w.writerow((index, cid, pid, activity))

    return drugs, proteins, interactions

from __future__ import print_function
import os, csv
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

    ### bidd.nus ###
    ### bidd.nus/drugDataset.txt ###
    # Files do have column headers (need to skip first line).
    # Three of the entries in drugDataset.txt don't have PubChemCIDs.
    #  pubchem_smiles_to_cid is able to find PubChemCIDs for all three of them.
    print('Processing bidd.nus/drugDataset.txt...')
    drugs = csv.reader(open('data/bidd.nus/drugDataset.txt', 'rb'), delimiter='\t')
    drugs.next() # Skip column header line.
    count = 0 # Keep track of how many entries successfully processed.
    for ttddrugid, drugname, pubchem_cid, smiles, inchi in drugs:
        assert ttddrugid not in ttd_drugs # Assert that there are no duplicate entries.
        if not pubchem_cid: # Look up PubChemCID based on SMILES if it's not there.
            pubchem_cid = pubchem_smiles_to_cid(smiles)
        assert pubchem_cid
        ttd_drugs[ttddrugid] = (drugname, pubchem_cid, smiles, inchi)
        count += 1
    print('Processed %d lines' % count)
    ### bidd.nus/targetDataset.txt ###
    # Many of the entries have residue numbers.
    # One entry does not have a UniProtID at all.
    #  Manual searches have not found a matching FASTA.
    # The sequences are listed without FASTA header.
    print('Processing bidd.nus/targetDataset.txt...')
    targets = csv.reader(open('data/bidd.nus/targetDataset.txt', 'rb'), delimiter='\t')
    targets.next() # Skip column header line.
    count = 0 # Keep track of how many entries successfully processed.
    for ttdtargetid, uniprotid, targetname, sequence in targets:
        assert ttdtargetid not in ttd_targets
        try:
            assert uniprotid
        except AssertionError as e:
            print('Doesn\'t have uniprotid', (ttdtargetid, uniprotid, targetname, sequence))
            continue
        assert sequence
        ttd_targets[ttdtargetid] = (uniprotid, targetname, sequence)
        count += 1
    print('Processed %d lines' % count)
    ### bidd.nus/interactionDataset.txt ###
    print('Processing bidd.nus/interactionDataset.txt...')
    interactions = csv.reader(open('data/bidd.nus/interactionDataset.txt', 'rb'), delimiter='\t')
    interactions.next() # Skip column header line.
    count = 0 # Keep track of how many entries successfully processed.
    for ttd_drugid, ttd_targetid, targetname, actions in interactions:
        assert ttd_drugid in ttd_drugs
        try:
            assert ttd_targetid in ttd_targets
        except AssertionError as e:
            print('Unknown target id', ttd_targetid, (ttd_drugid, ttd_targetid, targetname, actions))
            continue
        assert (ttd_drugid, ttd_targetid) not in ttd_interactions
        ttd_interactions[(ttd_drugid, ttd_targetid)] = (targetname, actions)
        count += 1
    print('Processed %d lines' % count)
    return ttd_drugs, ttd_targets, ttd_interactions

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
    drugs = csv.reader(open('data/drugbankversionmonth12Final/drugDataset.txt', 'rb'), delimiter='\t')
    count = 0 # Keep track of how many entries successfully processed.
    for db_drugid, drug_name, smiles, inchi, pubchem_cid in drugs:
        assert db_drugid not in db_drugs
        assert smiles
        if not pubchem_cid: # Look up PubChemCID based on SMILES if it's not there.
            #print('Missing pubchem_cid', db_drugid, drug_name, smiles, inchi, pubchem_cid)
            #pubchem_cid = pubchem_smiles_to_cid(smiles)
            continue
        db_drugs[db_drugid] = (drug_name, smiles, inchi, pubchem_cid)
        count += 1
    print('Processed %d lines' % count)

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
    interactions = csv.reader(open('data/drugbankversionmonth12Final/interactionsDataset.txt', 'rb'), delimiter='\t')
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


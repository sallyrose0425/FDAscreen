import pandas as pd
import json
import re
from collections import OrderedDict
import os
import sys
from time import sleep
import requests
# Replace all the urllib garbage with requests
if sys.version_info.major == 2:
    from urllib import urlopen, urlencode, quote
elif sys.version_info.major == 3:
    from urllib.parse import urlencode, quote
    from urllib.request import urlopen

from_types = ["ACC","ID","UPARC","NF90","NF100","GENENAME","EMBL_ID","EMBL","P_ENTREZGENEID","P_GI",
             "PIR","REFSEQ_NT_ID","P_REFSEQ_AC","UNIGENE_ID","PDB_ID","DISPROT_ID","BIOGRID_ID",
             "DIP_ID","MINT_ID","STRING_ID","CHEMBL_ID","DRUGBANK_ID","GUIDETOPHARMACOLOGY_ID",
             "SWISSLIPIDS_ID","ALLERGOME_ID","ESTHER_ID","MEROPS_ID","MYCOCLAP_ID","PEROXIBASE_ID",
             "REBASE_ID","TCDB_ID","BIOMUTA_ID","DMDM_ID","WORLD_2DPAGE_ID","DNASU_ID","ENSEMBL_ID",
             "ENSEMBL_PRO_ID","ENSEMBL_TRS_ID","ENSEMBLGENOME_ID","ENSEMBLGENOME_PRO_ID",
             "ENSEMBLGENOME_TRS_ID","GENEDB_ID","P_ENTREZGENEID","KEGG_ID","PATRIC_ID","UCSC_ID",
             "VECTORBASE_ID","WBPARASITE_ID","ARACHNOSERVER_ID","ARAPORT_ID","CCDS_ID","CGD",
             "CONOSERVER_ID","DICTYBASE_ID","ECHOBASE_ID","ECOGENE_ID","EUHCVDB_ID","EUPATHDB_ID",
             "FLYBASE_ID","GENECARDS_ID","GENEREVIEWS_ID","H_INVDB_ID","HGNC_ID","HPA_ID",
             "LEGIOLIST_ID","LEPROMA_ID","MAIZEGDB_ID","MGI_ID","MIM_ID","NEXTPROT_ID","ORPHANET_ID",
             "PHARMGKB_ID","POMBASE_ID","PSEUDOCAP_ID","RGD_ID","SGD_ID","TUBERCULIST_ID","WORMBASE_ID",
             "WORMBASE_PRO_ID","WORMBASE_TRS_ID","XENBASE_ID","ZFIN_ID","EGGNOG_ID","GENETREE_ID",
             "HOGENOM_ID","HOVERGEN_ID","KO_ID","OMA_ID","ORTHODB_ID","TREEFAM_ID","BIOCYC_ID",
             "REACTOME_ID","UNIPATHWAY_ID","CLEANEX_ID","COLLECTF_ID","CHITARS_ID","GENEWIKI_ID",
             "GENOMERNAI_ID"]

def check_uniprot_id(protein_id):
    url = 'http://www.uniprot.org/uniprot/%s.txt' % protein_id
    response = urlopen(url)
    if response.code == 200:
        return True
    else:
        return False

def strip_fasta(s):
    return ''.join(s.split('\n')[1:])

def get_fasta(uniprot_pid):
    url = 'http://www.uniprot.org/uniprot/%s.fasta' % uniprot_pid
    response = urlopen(url)
    if response.code == 200:
        return strip_fasta(response.read())
    else:
        return None

def smiles_to_cid(smiles):
    request_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/identity/smiles/JSON'
    callback_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/%s/cids/JSON'
    response = requests.post(request_url, [('smiles', smiles)])
    response = json.loads(response.text)
    while 'Waiting' in response:
        sleep(1)
        listkey = response['Waiting']['ListKey']
        raw_response = urlopen(callback_url % listkey).read()
        response = json.loads(raw_response)
    if not 'IdentifierList' in response:
        return None
    return response['IdentifierList']['CID'][0]

def name_to_cid(name):
    request_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/synonyms/JSON'
    response = requests.get(request_url % name)
    if response.status_code != 200:
        return None
    response_dict = json.loads(raw_response.text)
    return response_dict['InformationList']['Information'][0]['CID']

def cid_to_names(cid):
    request_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/JSON'
    response = requests.get(request_url % cid)
    if response.status_code != 200:
        return None
    response_dict = json.loads(response.text)
    return response_dict['InformationList']['Information'][0]['Synonym']


def get_uniprot_ids(protein_id):
    base_url = 'http://www.uniprot.org/mapping?%s'
    params = {'to': '',
              'format': 'tab',
              'query': protein_id}
    results = []
    for from_type in from_types:
        params['from'] = from_type
        url = base_url % urlencode(params)
        response = urlopen(url)
        content = response.read()
        lines = content.split(b'\n')
        if lines[1]:
            from_id, to_id = lines[1].split(b'\t')
            print(from_type,to_id)
            results.append((from_type, to_id))
    return results

def get_pubchem_json(protein_id):
    base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/protein/%s/JSON'

def query_cid_interactions(cid):
    base_url = 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?%s'
    query = [OrderedDict([('download', ['activity', 'acvalue', 'acname','targetname',
                                        'targeturl','aidname', 'aid', 'sid', 'cid']),
                          ('collection',   'bioactivity'),
                          ('where',        {'ands': [{'cid': str(cid)}] }),
                          ('order',        ['acvalue,asc']),
                          ('start',        1),
                          ('limit',        1000000),
                          ('nullatbottom', 1)
                         ]),
             OrderedDict([('histogram', 'activity'), ('bincount', 10000)])]
    params = {'infmt': 'json', 'outfmt': 'json', 'query': json.dumps(query)}
    url = base_url % urlencode(params)
    table = pd.read_csv(url, na_filter=False, engine='python')
    return table

def pubchem_to_uniprot(pubchem_pid):
    pid = fasta = None
    pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/protein/%s' % pubchem_pid
    pubchem_html = str(urlopen(pubchem_url).read())
    matches = re.findall('<meta name="pubchem_uid_value" content="([^"]+)">', pubchem_html)
    if matches:
        pid = matches[0]
        fasta = get_fasta(pid)
    if not pid or not fasta:
        print("Couldn't get pid from %s!" % pubchem_pid)
    return pid, fasta

def main():
    interactions = {'cid':[], 'pid':[], 'activity':[]}
    compounds = {'cid':[], 'smiles':[]}
    proteins = {'pid':[], 'fasta':[]}

    drugs = pd.read_csv('data/drugDataset.txt', delimiter='\t',
                         names=['dbid', 'name', 'smiles', 'inchi', 'cid'],
                         na_filter=False)

    for index, cid, smiles in drugs[['cid','smiles']].itertuples():
        if not cid and smiles:
            cid = smiles_to_cid(smiles)
        if not cid:
            print("Couldn't get cid!")
            continue
        print(cid)
        if not smiles:
            print("Doesn't have smiles!")
            continue
        if cid not in compounds['cid']:
            compounds['cid'].append(cid)
            compounds['smiles'].append(smiles)
        if cid in interactions['cid']:
            print('Skipping cid', cid)
            continue
        table = query_cid_interactions(cid)
        if 'activity' not in table or 'targeturl' not in table:
            print('Unable to get interactions!')
            continue
        for index, activity, target_url in table[['activity', 'targeturl']].itertuples():
            if activity in ('Active', 'Inactive') and '/protein/' in target_url:
                pubchem_pid = target_url.split('/')[-1]
                pid, fasta = pubchem_to_uniprot(pubchem_pid)
                if not pid or not fasta:
                    continue
                if pid not in proteins['pid']:
                    proteins['pid'].append(pid)
                    proteins['fasta'].append(fasta)
                interactions['cid'].append(cid)
                interactions['pid'].append(pid)
                interactions['activity'].append(int(activity=='Active'))
                print(cid, pid, activity=='Active')
    return compounds, proteins, interactions

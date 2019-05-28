import pandas as pd
import urllib
import json
import re
from collections import OrderedDict
import os

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
    try:
        response = urllib.request.urlopen(url)
        return True
    except Exception:
        return False

def get_uniprot_ids(protein_id):
    base_url = 'http://www.uniprot.org/mapping?%s'
    params = {'to': '',
              'format': 'tab',
              'query': protein_id}
    results = []
    for from_type in from_types:
        params['from'] = from_type
        url = base_url % urllib.parse.urlencode(params)
        response = urllib.request.urlopen(url)
        content = response.read()
        lines = content.split(b'\n')
        if lines[1]:
            from_id, to_id = lines[1].split(b'\t')
            print(from_type,to_id)
            results.append((from_type, to_id))
    return results

def get_pubchem_json(protein_id):
    base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/protein/%s/JSON'

def get_cid_interactions(cid):
    if os.path.exists('tmp/skip_pids.txt'):
        skip_pids = eval(open('tmp/skip_pids.txt', 'r').read())
    else:
        skip_pids = []
    if os.path.exists('tmp/pid_map.txt'):
        pid_map = eval(open('tmp/pid_map.txt', 'r').read())
    else:
        pid_map = {}
    seen = []
    interactions = []
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
    url = base_url % urllib.parse.urlencode(params)
    table = pd.read_csv(url, na_filter=False)
    if 'activity' not in table or 'targeturl' not in table:
        return []
    ids = []
    for i, (activity, target_url) in table[['activity', 'targeturl']].iterrows():
        if activity in ('Active', 'Inactive') and '/protein/' in target_url:
            pid = target_url.split('/')[-1]
            if pid in skip_pids:
                continue
            if pid not in pid_map:
                pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/protein/%s' % pid
                pubchem_html = str(urllib.request.urlopen(pubchem_url).read())
                matches = re.findall('<meta name="pubchem_uid_value" content="([^"]+)">', pubchem_html)
                if matches:
                    if check_uniprot_id(matches[0]):
                        pid_map[pid] = matches[0]
                        pid = matches[0]
                    else:
                        skip_pids.append(matches[0])
                        continue
                else:
                    continue
            else:
                pid = pid_map[pid]
            if (cid,pid) not in seen:
                seen.append((cid,pid))
                interactions.append((cid, pid, int(activity=='Active')))
    with open('tmp/skip_pids.txt', 'w') as f:
        f.write(repr(skip_pids))
    with open('tmp/pid_map.txt', 'w') as f:
        f.write(repr(pid_map))
    return interactions


def main():
    if not os.path.exists('tmp'):
        os.mkdir('tmp')
    if os.path.exists('tmp/skip_cids.txt'):
        skip_cids = eval(open('tmp/skip_cids.txt', 'r').read())
    else:
        skip_cids = []
    if os.path.exists('data/interactions.csv'):
        df = pd.read_csv('data/interactions.csv', index_col=0, na_filter=False,
                         dtype={'cid':str,'pid':str,'activity':int})
    else:
        df = pd.DataFrame(columns=['cid', 'pid', 'activity'])

    drugs = pd.read_csv('data/drugDataset.txt', delimiter='\t', 
                         names=['drugbankId', 'drugName', 'drugSMILES', 'drugInChI', 'drugPubChemCompound'], 
                         na_filter=False)

    for cid in drugs['drugPubChemCompound']:
        if cid:
            if cid in skip_cids:
                print('Skipping cid', cid)
            elif cid in df['cid'].values:
                print('Skipping cid', cid)
                skip_cids.append(cid)
            else:
                print(cid)
                interactions = get_cid_interactions(cid)
                print(interactions)
                rows = pd.DataFrame(interactions, columns=['cid', 'pid', 'activity'], dtype=object)
                df = df.append(rows, ignore_index=True)
                skip_cids.append(cid)
                with open('tmp/skip_cids.txt', 'w') as f:
                    f.write(repr(skip_cids))
                df.to_csv('data/interactions.csv')

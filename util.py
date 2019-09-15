import json, requests, csv, re
import pandas as pd
from time import sleep
from dl import query_cid_interactions
# Replace all the urllib garbage with requests
import sys
if sys.version_info.major == 2:
    from urllib import urlopen, urlencode, quote
elif sys.version_info.major == 3:
    from urllib.parse import urlencode, quote
    from urllib.request import urlopen

# Come up with a better name for this function.
# Reads the drugbankversionmonth12Final/proteinDataset.txt file.
def read_broken_file(filename):
    with open(filename) as f:
        records = []
        record = line = f.readline()
        while line:
            line = f.readline()
            if '\t' in line or not line:
                records.append(record.split('\t'))
                record = ''
            record += line
    return records

def pubchem_get_canonical_pid(pid):
    s = requests.Session()
    a = requests.adapters.HTTPAdapter(max_retries=5)
    s.mount('https://', a)
    response = s.get('https://pubchem.ncbi.nlm.nih.gov/protein/%s' % pid)
    html = response.text
    pid1_matches = re.findall('<meta name="pubchem_uid_value" content="([^\"]*)">', html)
    pid2_matches = re.findall('<link rel="canonical" href="https://pubchem.ncbi.nlm.nih.gov/protein/([^\"]*)"/>', html)
    assert pid1_matches and pid2_matches
    pid1, pid2 = pid1_matches[0], pid2_matches[0]
    assert pid1 and pid2 and pid1 == pid2
    return pid1

# Assumes pid is a canonical pubchem protein id.
def pubchem_pid_to_fasta(pid):
    data_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/protein/%s/JSON'
    s = requests.Session()
    a = requests.adapters.HTTPAdapter(max_retries=5)
    s.mount('https://', a)
    response = s.get(data_url % pid)
    data = json.loads(response.text)
    fasta = ''
    for section in data['Record']['Section']:
        if section['TOCHeading'] == 'Sequence':
            assert not fasta, 'Found multiple fasta sections'
            fasta = section['Information'][0]['Value']['StringWithMarkup'][1]['String']
    assert fasta, "Didn't find a fasta section"
    return fasta

# Strengthen error checking.
# Compare to searches based on drug name.
def pubchem_smiles_to_cid(smiles):
    try:
        with open('tmp/pubchem_smiles_to_cid.cache', 'r') as f:
            cache = eval(f.read())
    except FileNotFoundError as e:
        cache = {}
    if smiles in cache:
        return cache[smiles]
    request_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/identity/smiles/JSON'
    callback_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/%s/cids/JSON'
    response = requests.post(request_url, [('smiles', smiles)])
    response = json.loads(response.text)
    while 'Waiting' in response:
        sleep(1)
        listkey = response['Waiting']['ListKey']
        raw_response = urlopen(callback_url % listkey).read()
        response = json.loads(raw_response)
    assert 'IdentifierList' in response
    cid = response['IdentifierList']['CID'][0]
    cache[smiles] = cid
    with open('tmp/pubchem_smiles_to_cid.cache', 'w') as f:
        f.write(repr(cache))
    return cid

def strip_fasta(s):
    return ''.join(s.split('\n')[1:])

import pandas as pd
import importlib
import omnipath as op
from tqdm import tqdm
import re
import id_conversion

hgnc_symbols = pd.read_csv("././data/interim/HGNC symbols.txt", sep=',', header=0)
omnipath = op.interactions.AllInteractions.get(
    directed=False, organism='human')

complex_a = list(omnipath.loc[omnipath['source'].str.startswith(
    'COMPLEX:', na=False)]['source'])
complex_b = list(omnipath.loc[omnipath['target'].str.startswith(
    'COMPLEX:', na=False)]['target'])
complex_ids = set(complex_a + complex_b)

omnipath_dict = omnipath.to_dict('records')
row_list = []
columns = ['source', 'target']
for row in omnipath_dict:
    for column in columns:
        comp = re.findall(r'(?<=[:_])([a-zA-Z0-9]*)', row[column])
        if len(comp) > 0:
            row[column] = comp
        else:
            row[column] = [row[column]]  # list to iterate in the next step

    for prot1 in row[columns[0]]:
        for prot2 in row[columns[1]]:
            row_list.append(
                {columns[0]: prot1, columns[1]: prot2, 'is_directed': row['is_directed']})
omnipath = pd.DataFrame(row_list).drop_duplicates()

omnipathIDs = set(list(omnipath['source'])+list(omnipath['target']))

for id in tqdm(omnipathIDs):
    hgnc_id = hgnc_symbols[hgnc_symbols['UniProt ID(supplied by UniProt)']
                           == id]['Approved symbol'].values
    if len(hgnc_id) > 0:
        omnipath.loc[omnipath['source'] == id, 'HGNC A'] = hgnc_id[0]
        omnipath.loc[omnipath['target'] == id, 'HGNC B'] = hgnc_id[0]

omnipath_missing_ids_a = list(
    omnipath[omnipath['HGNC A'].isnull()]['source'].unique())
omnipath_missing_ids_b = list(
    omnipath[omnipath['HGNC B'].isnull()]['target'].unique())
omnipath_missing_ids = set(omnipath_missing_ids_a + omnipath_missing_ids_b)

mimat = []
non_mimat_missing_ids = []
for id in omnipath_missing_ids:
    if id.startswith('MIMAT'):
        mimat.append(id)
    else:
        non_mimat_missing_ids.append(id)

uniprot_search_matches, uniprot_search_misses = id_conversion.uniprot_search(
    non_mimat_missing_ids, 'ACC+ID', 'GENENAME')
uniprot_confirmed_hgnc, uniprot_ids_no_match = id_conversion.confirm_hgnc(
    list(uniprot_search_matches.keys()))
non_uniprot_matches, no_match = id_conversion.confirm_hgnc(
    uniprot_search_misses)

uniprot_confirmed_hgnc_df = pd.DataFrame.from_dict(
    uniprot_confirmed_hgnc, orient='index')
non_uniprot_matches_df = pd.DataFrame.from_dict(
    non_uniprot_matches, orient='index')
hgnc_symbols_df = pd.concat(
    [uniprot_confirmed_hgnc_df, non_uniprot_matches_df])

for id in tqdm(non_mimat_missing_ids):
    hgnc_id = hgnc_symbols_df[hgnc_symbols_df.isin([id]).any(axis=1)].index
    if len(hgnc_id) > 0:
        omnipath.loc[omnipath['source'] == id, 'HGNC A'] = hgnc_id[0]
        omnipath.loc[omnipath['target'] == id, 'HGNC B'] = hgnc_id[0]

mirbase_matches = id_conversion.scrape_mirbase(mimat)
mir_matches_df = pd.DataFrame.from_dict(mirbase_matches, orient='index')
for id in mimat:
    hgnc_id = mir_matches_df[mir_matches_df.isin([id]).any(axis=1)].index
    if len(hgnc_id) > 0:
        omnipath.loc[omnipath['source'] == id, 'HGNC A'] = hgnc_id[0]
        omnipath.loc[omnipath['target'] == id, 'HGNC B'] = hgnc_id[0]

omnipath = omnipath[~((omnipath['HGNC A'].isnull()) |
                      (omnipath['HGNC B'].isnull()))]
omnipath.to_csv('././data/interim/omnipath.csv', header=-1, index=False)

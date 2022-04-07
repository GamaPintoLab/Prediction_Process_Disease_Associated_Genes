import pandas as pd
import numpy as np

apid = pd.read_csv("././data/interim/apid.csv", sep=',',
                   header=0, names=['HGNC_A', 'HGNC_B'])
apid['apid'] = True
apid.dropna(inplace=True)

dorothea = pd.read_csv("././data/interim/dorothea.csv",
                       sep=',', header=0, names=['HGNC_A', 'HGNC_B'])
dorothea['dorothea'] = True

omnipath = pd.read_csv("././data/interim/omnipath.csv", sep=',', header=0,
                       names=['source', 'target', 'is_directed', 'HGNC_A', 'HGNC_B'])
omnipath['omnipath'] = True

huri = pd.read_csv("././data/interim/HuriHGNC_IDs.csv", sep=',', header=0)
huri['huri'] = True

ppi_interactions = pd.merge(apid[['HGNC_A', 'HGNC_B', 'apid']], dorothea[[
                            'HGNC_A', 'HGNC_B', 'dorothea']], on=['HGNC_A', 'HGNC_B'], how='outer')
ppi_interactions = pd.merge(ppi_interactions, huri[['HGNC_A', 'HGNC_B', 'huri']], on=[
                            'HGNC_A', 'HGNC_B'], how='outer')
ppi_interactions = pd.merge(ppi_interactions, omnipath[[
                            'HGNC_A', 'HGNC_B', 'omnipath', 'is_directed']], on=['HGNC_A', 'HGNC_B'], how='outer')
ppi_interactions.fillna(False, inplace=True)
ppi_interactions.to_csv(
    '././data/processed/ppis/ppi_interactions.csv', header=-1, index=False)

ppi_interactions = pd.read_csv('././data/processed/ppis/ppi_interactions.csv')
apid_huri_ppis = ppi_interactions[(ppi_interactions['apid'] == True) | (
    ppi_interactions['huri'] == True)]
apid_huri_ppis.to_csv(
    '././data/processed/ppis/apid_huri_ppis.csv', index=False, header=-1)

disgenet = pd.read_csv(
    "././data/raw/curated_gene_disease_associations.tsv", sep='\t', header=0)
disgenet.to_csv('././data/interim/disgenet.csv', header=-1, index=False)

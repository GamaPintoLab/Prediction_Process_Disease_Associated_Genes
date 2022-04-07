import requests
import re
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.common.exceptions import NoSuchElementException, TimeoutException
from tqdm.notebook import tqdm


def uniprot_search(ids, original_id, final_id, id_range=1000):
    # Uses UniProt API in order to convert a set of ids to HGNC GeneName (or other).
    #
    # INPUT:
    #   -set of ids to convert
    #   -id format that the proteins are in
    #   -id format that we want to convert to
    #   -number of ids per request
    #
    # RETURNS: one dict with the final ids as keys and the original ones as values and a list with the ids that missed

    id_map = {}
    no_match = []
    ids = list(ids)
    for pos in tqdm(range(0, len(ids), id_range)):
        idset = ids[pos:pos+id_range]
        id = ' '.join(idset)
        params = {
            'from': original_id,
            'to': final_id,
            'format': 'tab',
            'query': id}
        response = requests.get(
            'https://www.uniprot.org/uploadlists/',
            params=params)
        resulting_codes = list(filter(None, re.split(
            '\t|\n|To|From', str(response.text))))
        for pos in range(0, len(resulting_codes), 2):
            if resulting_codes[pos] in idset:
                idset.remove(resulting_codes[pos])
            if resulting_codes[pos+1] not in id_map:
                id_map[resulting_codes[pos+1]] = [resulting_codes[pos]]
            else:
                id_map[resulting_codes[pos+1]].append(resulting_codes[pos])
        no_match.extend(idset)
    return id_map, no_match


def genecards_search(proteins, idtype='uniprot'):
    # Uses HGNC's REST web-service in order to convert uniprot, ensembl, ncbi or old HGNC IDs to current HGNC IDs.
    #
    # INPUT:
    #   -set of ids to convert
    #   -id format that the proteins are in
    #
    # RETURNS: one dict with the HGNC ids as keys and the original ones as values and a list with the ids that missed

    idtypes = {'uniprot': 'uniprot_ids', 'ensembl': 'ensembl_gene_id',
               'ncbi': 'entrez_id', 'alias symbol': 'alias_symbol'}
    id_map = {}
    no_match = []
    for id in tqdm(proteins):
        url = "http://rest.genenames.org/search/"+idtypes[idtype]+"/"+id
        response = requests.get(url, headers={'Accept': 'application/json'})
        json_info = response.json()
        if len(json_info['response']['docs']) == 1:
            for i in json_info['response']['docs']:
                hgnc = i['symbol']
                print(id, hgnc)
                if hgnc not in id_map.keys():
                    id_map[hgnc] = [id]
                else:
                    id_map[hgnc].append([id])
        else:
            print(id)
            no_match.append(id)
    return id_map, no_match


def ncbi_search(protein_ids):
    # Searches NCBI web-service to convert Ensembl and NCBI IDs to their HGNC GeneName
    #
    # INPUT:
    #   -set of ids to convert
    #
    # RETURNS: one dict with the HGNC ids as keys and the original ones as values
    id_map = {}
    for id in tqdm(protein_ids):
        response = requests.get(
            "https://www.ncbi.nlm.nih.gov/gene/?term=" + id + "&report=full_report&format=text")
        if len(response.text) <= 172:
            print(id)
        else:
            data = response.text.splitlines()
            for line in data:
                if line.startswith('Official Symbol:') and line.endswith('(provided by HGNC)'):
                    hgnc_code = line.split(' ')[2]
                    id_map[hgnc_code] = [id]
                    print(id, hgnc_code)
                    continue
            if [id] not in id_map.values():
                print(id)
    return id_map


def confirm_hgnc(protein_list):
    # Uses HGNC's REST web-service in order to confirm is the HGNC IDs obtained in other sources are the current ones,
    # or if they're aliases, older symbols or names of the gene while retriving their updated ID.
    #
    # INPUT:
    #   -set of ids to convert
    #
    # RETURNS: one dict with the HGNC ids as keys with the original ones as values, and a list with the unresolved ones.
    idtypes = {'symbol': 'symbol', 'previous symbol': 'prev_symbol', 'name': 'name', 'alias symbol': 'alias_symbol',
               'previous name': 'prev_name',  'alias names': 'alias_name'}
    resolved_ids = {}
    for idtype in idtypes.values():
        if len(protein_list) == 0:
            break
        print('Getting results for {}'.format(idtype))
        print('Testing {} proteins'.format(len(protein_list)))
        print('Number of resolved IDs: {}'.format(len(resolved_ids)))
        proteins_to_remove = []
        for protein in tqdm(protein_list):
            if '/' in protein:

                if protein not in resolved_ids.keys():
                    resolved_ids[protein] = [protein]
                else:
                    resolved_ids[protein].extend([protein])
                proteins_to_remove.append(protein)
            else:
                url = "http://rest.genenames.org/search/"+idtype+"/"+protein
                response = requests.get(
                    url, headers={'Accept': 'application/json'})
                json_info = response.json()
                if len(json_info['response']['docs']) == 1 and idtype == 'symbol':
                    for i in json_info['response']['docs']:
                        hgnc = i['symbol']
                        print(protein, hgnc)
                        if hgnc not in resolved_ids.keys():
                            resolved_ids[hgnc] = [protein]
                        else:
                            resolved_ids[hgnc].extend([protein])
                        proteins_to_remove.append(protein)
                elif 0 < len(json_info['response']['docs']) < 4 and idtype != 'symbol':
                    protein_info = json_info['response']['docs'][0]
                    hgnc = protein_info['symbol']
                    print(protein, hgnc)
                    if hgnc not in resolved_ids.keys():
                        resolved_ids[hgnc] = [protein]
                    else:
                        resolved_ids[hgnc].extend([protein])
                    proteins_to_remove.append(protein)
        protein_list = [
            protein for protein in protein_list if protein not in proteins_to_remove]
    return resolved_ids, protein_list


def scrape_mirbase(proteins):
    # Uses selenium to scrape the miRBase in order to convert miRBase Accession Number IDs to HGNC Gene Names.
    #
    # INPUT:
    #   -set of ids to convert
    #
    # RETURNS: one dict with the HGNC ids as keys with the original ones as values.
    matches = {}
    options = Options()
    options.headless = True
    driver = webdriver.Firefox(options=options)
    for id in tqdm(proteins):
        url = "https://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=" + id
        driver.get(url)
        try:
            element = driver.find_element_by_xpath(
                '//*[@id="symbol"]/a')
            element = element.text
        except (NoSuchElementException, TimeoutException):
            print("\n\nNO MATCH FOUND IN {}\n\n".format(id))
        finally:
            if element.startswith('HGNC'):
                hgnc_id = element.split(':')[1]
                print(id, hgnc_id)
                if hgnc_id in matches:
                    matches[hgnc_id].append(id)
                else:
                    matches[hgnc_id] = [id]
            else:
                print("\n\nNO MATCH FOUND IN {}\n\n".format(id))
    driver.quit()
    return matches

from multiprocessing import process
from os import close
import pandas as pd
import numpy as np
import math
from tqdm.notebook import tqdm
from scipy.stats import hypergeom
from collections import Counter
from itertools import combinations
from igraph import Graph
from itertools import chain
from ast import literal_eval


def hypergeometric_test(graph, proteins_by_process_df, adjacency_matrix):
    # Computation of the negative log10 p-value given by an hypergeometric test.
    #
    # INPUT:
    #   - graph with the PPI network
    #   - dataframe with proteins of each process/disease
    #   - adjacency matrix of the graph
    #
    # RETURNS: dataframe with -log10(p-value) values for every protein in every process/disease.
    proteins_neglogp_byprocess = {}
    protein_indices = graph.vs.indices
    protein_names = graph.vs['name']
    n_proteins = len(protein_names)
    proteins_by_process_dict = proteins_by_process_df.to_dict('index')
    for index in tqdm(protein_indices):
        protein = protein_names[index]
        processes_neglogp = protein_hypertests(
            index, proteins_by_process_dict, adjacency_matrix, n_proteins)
        proteins_neglogp_byprocess[protein] = processes_neglogp
    return proteins_neglogp_byprocess


def protein_hypertests(protein_index, proteins_by_process_dict, adjacency_matrix, n_proteins):
    # Computation of the negative log10 p-value given by an hypergeometric test for a given protein in all processes/diseases.
    #
    # INPUT:
    #   - index of the query protein
    #   - dict with proteins present in each process/disease
    #   - adjacency matrix of the graph
    #   - number of proteins in the network
    #
    # RETURNS: dict with -log10(p-value) values for the target protein in every process/disease.
    protein_process_neglogpvalue = {}
    protein_row = adjacency_matrix[protein_index]
    n_ppi_target_protein = sum(protein_row)
    for index, values in proteins_by_process_dict.items():
        num_neighbors = 0
        process_protein_indexes = list(values['protein_index'])
        process_n_proteins = values['n_proteins']
        for process_protein_index in process_protein_indexes:
            if process_protein_index > adjacency_matrix.shape[0]:
                continue
            if protein_row[process_protein_index] == 1:
                num_neighbors += 1
        neg_logp = -math.log(hypergeom.pmf(num_neighbors,
                             n_proteins-1, process_n_proteins, n_ppi_target_protein))
        protein_process_neglogpvalue[values['process']
                                     ] = neg_logp
    return protein_process_neglogpvalue


def shortest_paths(graph, processes_df):
    # Computation of the shortest paths betweeen every protein and a process/disease protein.
    #
    # INPUT:
    #   - graph with the PPI network
    #   - dataframe with processes/diseases and their respective proteins.
    #
    # RETURNS: dataframe with shortest paths between any protein (rows) and a process/disease proteins (columns).
    proteins_shortest_paths = {}
    all_process_proteins = list(set(
        [protein for protein_list in processes_df['proteins_ids'].values for protein in protein_list if protein in graph.vs['name']]))
    target_proteins = graph.vs['name']
    for protein in all_process_proteins:
        if protein not in target_proteins:
            all_process_proteins.remove(protein)
    shortest_paths = graph.shortest_paths(
        source=all_process_proteins, target=target_proteins, mode='all')
    for source_protein_index in tqdm(range(len(all_process_proteins))):
        source_protein = all_process_proteins[source_protein_index]
        proteins_shortest_paths[source_protein] = {}
        for target_protein_index in range(len(target_proteins)):
            target_protein = target_proteins[target_protein_index]
            short_path_len = shortest_paths[source_protein_index][target_protein_index]
            proteins_shortest_paths[source_protein][target_protein] = short_path_len
    process_shortest_paths_df = pd.DataFrame.from_dict(proteins_shortest_paths)
    process_shortest_paths_df = process_shortest_paths_df.rename(
        index=dict(zip(list(process_shortest_paths_df.index), graph.vs['name'])))
    return process_shortest_paths_df


def closeness(shortest_paths_df, process_df):
    # Computation of the closeness of every protein to the proteins of each process/disease.
    #
    # INPUT:
    #   - dataframe of shortest paths
    #   - dataframe with processes/diseases and their respective proteins.
    #
    # RETURNS: dataframe with closeness of a protein in each process/disease.
    shortest_paths = shortest_paths_df.to_dict('dict')
    process_proteins_closeness = {}
    for index, process in tqdm(process_df.iterrows(), total=process_df.shape[0]):
        n_process_proteins = len(process['proteins_ids'])
        process_id = process['process']
        process_proteins = {
            proteins: shortest_paths[proteins] for proteins in process['proteins_ids'] if proteins in shortest_paths_df.index}
        process_proteins_values = process_proteins.values()
        shortest_sum = Counter()
        for protein in process_proteins_values:
            shortest_sum.update(protein)
        shortest_sum = dict(shortest_sum)
        #closeness = {protein: short_sum / n_process_proteins for protein, short_sum in shortest_sum.items()}
        closeness = {protein: n_process_proteins /
                     short_sum for protein, short_sum in shortest_sum.items()}
        process_proteins_closeness[process_id] = closeness
    return process_proteins_closeness


def betweenness(shortest_paths, process_df, graph):
    # Computation of the binary betweenness of every protein to the proteins of each process/disease.
    # For every shortest paths between two proteins of the same process, checks what are the middle proteins (connecting both proteins) and keeps track of the amount of times they connect to process proteins.
    #
    # INPUT:
    #   - dataframe of shortest paths
    #   - dataframe with processes/diseases and their respective proteins.
    #   - graph with PPIs
    #
    # RETURNS: dataframe with betweenness of a protein in each process/disease.
    process_proteins_betweenness = {}
    for index, process in tqdm(process_df.iterrows(), total=process_df.shape[0]):
        process_proteins = [
            protein for protein in process['proteins_ids'] if protein in shortest_paths.index]
        process_id = process['process']
        process_proteins_pairs = combinations(process_proteins, 2)
        n_possible_paths = int(
            (len(process_proteins)*(len(process_proteins)-1))/2)
        process_proteins_betweenness[process_id] = {}
        shortest_sum = Counter()
        for pair in process_proteins_pairs:
            pair_shortest_path = shortest_paths.at[pair[0], pair[1]]
            intermediate_paths = shortest_paths.loc[:, [
                pair[0], pair[1]]].sum(axis=1, skipna=False)
            middle_proteins = list(
                intermediate_paths[intermediate_paths == pair_shortest_path].index)
            shortest_sum.update(middle_proteins)
        process_proteins_betweenness[process_id] = {protein: short_sum /
                                                    n_possible_paths for protein, short_sum in shortest_sum.items()}
        for protein in graph.vs['name']:
            if protein not in process_proteins_betweenness[process_id].keys():
                process_proteins_betweenness[process_id][protein] = 0
    return process_proteins_betweenness


def fraction_betweenness(processes_df, graph):
    # Computation of the binary betweenness of every protein to the proteins of each process/disease.
    # For every shortest paths between two proteins of the same process, checks what are the middle proteins (connecting both proteins) and keeps track of the amount of times they connect to process proteins.
    #
    # INPUT:
    #   - dataframe of shortest paths
    #   - dataframe with processes/diseases and their respective proteins.
    #   - graph with PPIs
    #
    # RETURNS: dataframe with betweenness of a protein in each process/disease.
    try:
        processes_df['protein_index'] = processes_df['protein_index'].apply(
            literal_eval)
        processes_df['proteins_ids'] = processes_df['proteins_ids'].apply(
            literal_eval)
    except ValueError:
        processes_df['protein_index'] = processes_df['protein_index']
        processes_df['proteins_ids'] = processes_df['proteins_ids']
    #processes_df = processes_df.iloc[:3,:]
    process_proteins = list(set(
        [protein for protein_list in processes_df['proteins_ids'].values for protein in protein_list if protein in graph.vs['name']]))
    n_combinations = int((len(process_proteins)*(len(process_proteins)-1))/2)
    processes_by_proteins = {}
    for k, v in processes_df[['process', 'proteins_ids']].set_index('process').to_dict('dict')['proteins_ids'].items():
        for x in v:
            if x in graph.vs['name']:
                processes_by_proteins.setdefault(
                    int(graph.vs.find(name=x)['id']), []).append(k)

    process_proteins_pairs = combinations(process_proteins, 2)
    protein_pairs = {}
    for pair in process_proteins_pairs:
        if pair[0] in protein_pairs:
            protein_pairs[pair[0]].append(pair[1])
        else:
            protein_pairs[pair[0]] = [pair[1]]

    pairs_shortest_paths = {}

    for protein, protein_list in tqdm(protein_pairs.items()):
        shortest_paths = graph.get_all_shortest_paths(protein, protein_list)
        for shortest_path in shortest_paths:

            if (shortest_path[0], shortest_path[-1]) not in pairs_shortest_paths:
                pairs_shortest_paths[(
                    shortest_path[0], shortest_path[-1])] = [shortest_path[1:-1]]
            else:
                pairs_shortest_paths[(
                    shortest_path[0], shortest_path[-1])].append(shortest_path[1:-1])
    process_proteins_betweenness = {}
    for process in processes_df['process']:
        process_proteins_betweenness[process] = {}

    for pair, protein_shortest_paths in tqdm(pairs_shortest_paths.items(), total=n_combinations):
        protein_a = pair[0]
        protein_b = pair[1]
        shortest_path_count = Counter(
            chain.from_iterable(protein_shortest_paths))
        shortest_path_count = {
            protein: n/len(protein_shortest_paths) for protein, n in shortest_path_count.items()}
        processes_in_pair = list(set(processes_by_proteins[protein_a]) & set(
            processes_by_proteins[protein_b]))
        for process in processes_in_pair:
            process_proteins_betweenness[process] = dict(
                Counter(process_proteins_betweenness[process])+Counter(shortest_path_count))

    for protein in range(len(graph.vs['name'])):
        for process in process_proteins_betweenness.keys():
            if protein not in process_proteins_betweenness[process]:
                process_proteins_betweenness[process][protein] = 0

    id_dict = {}
    for i in range(len(graph.vs['name'])):
        id_dict[i] = graph.vs['name'][i]

    betweenness = pd.DataFrame.from_dict(
        process_proteins_betweenness, orient='columns')
    betweenness['index'] = betweenness.index.to_series().map(id_dict)
    betweenness.set_index('index', inplace=True)
    betweenness.sort_index(inplace=True)
    return betweenness


def random_walk_restart(graph, process_df):
    # Computation of the random walks with restart of every protein to the proteins of each process/disease with igraph's pagerank algorithm.
    #
    # INPUT:
    #   - dataframe of shortest paths
    #   - dataframe with processes/diseases and their respective proteins.
    #
    # RETURNS: dataframe with random walks with restart scores of every protein in each process/disease.
    process_proteins_rwr = {}
    for index, process in tqdm(process_df.iterrows(), total=process_df.shape[0]):
        process_proteins = [
            protein for protein in process['proteins_ids'] if protein in graph.vs['name']]
        process_id = process['process']
        rwr_values = graph.personalized_pagerank(
            reset_vertices=process_proteins)
        process_proteins_rwr[process_id] = rwr_values
    return process_proteins_rwr


def multiple_metrics(ppis, process_df):
    # Allows for the computation of the hypergeometric test, shortest paths, closeness, betweenness, and random walks with restart in the reduced networks.
    #
    # INPUT:
    #   - numpy array with PPIs
    #   - dataframe with processes/diseases and their respective proteins.
    #
    # RETURNS: collection of dataframes with hypergeometric, closeness, betweenneess and random walk with restart scores.
    hyper_scores = []
    closeness_scores = []
    betweenness_scores = []
    rwr_scores = []
    fraction_betweenness_scores = []
    df_keys = []
    keys = 0
    for array in tqdm(ppis, total=len(ppis)):
        df_keys.append(keys)
        keys += 1
        print('Generating Graph')
        ppi_df = pd.DataFrame(array, columns=[
                              'HGNC_A', 'HGNC_B', 'apid', 'dorothea', 'huri', 'omnipath', 'is_directed'])
        ppi_df.replace(0, np.nan, inplace=True)
        ppi_df.dropna(axis=0, how='all')
        reduced_graph = Graph.DataFrame(ppi_df, directed=False)
        
        graph = reduced_graph.simplify()
        if not graph.is_connected():
            cluster = graph.clusters()
            graph = graph.induced_subgraph(cluster[0])
        graph.write_gml('../data/graph_apid_huri_80')
        graph = Graph.Read_GML("../data/graph_apid_huri_80")
        #print('Generating Confusion Matrix')
        #adj_matrix = graph.get_adjacency()
        #adj_matrix = np.array(adj_matrix.data)
        #print('Hypergeometric Test Computation:')
        #hyper_dict = hypergeometric_test(graph, process_df, adj_matrix)
        #hyper_df = pd.DataFrame.from_dict(hyper_dict).transpose()
        # hyper_scores.append(hyper_df)
        #print('Shortest Paths Computation:')
        #shortest_paths_df = shortest_paths(graph, process_df)
        #print('Closeness Computation:')
        #closeness_dict = closeness(shortest_paths_df, process_df)
        #closeness_df = pd.DataFrame.from_dict(closeness_dict)
        # closeness_scores.append(closeness_df)
        #print('Betweenness Computation:')
        #betweenness_dict = betweenness(shortest_paths_df, process_df, graph)
        #betweenness_df = pd.DataFrame.from_dict(betweenness_dict)
        #betweenness_df.fillna(value=0, inplace=True)
        # betweenness_scores.append(betweenness_df)
        print('Fraction Betweenness Computation:')
        fraction_betweenness_df = fraction_betweenness(process_df, graph)
        fraction_betweenness_df.fillna(value=0, inplace=True)
        fraction_betweenness_scores.append(fraction_betweenness_df)
        #print('Random-Walks with Restart Computation:')
        #rwr_dict = random_walk_restart(graph, process_df)
        #rwr_df = pd.DataFrame.from_dict(rwr_dict)
        # rwr_df.rename(index=dict(
        #    zip(list(rwr_df.index), graph.vs['name'])), inplace=True)
        # rwr_scores.append(rwr_df)
    # hyper_scores_df = pd.concat(
    #    hyper_scores, keys=df_keys, axis=0).reset_index(level=1)
    # closeness_df = pd.concat(
    #    closeness_scores, keys=df_keys, axis=0).reset_index(level=1)
    # betweenness_df = pd.concat(
    #    betweenness_scores, keys=df_keys, axis=0).reset_index(level=1)
    fraction_betweenness_df = pd.concat(
        fraction_betweenness_scores, keys=df_keys, axis=0).reset_index(level=1)
    #rwr_df = pd.concat(rwr_scores, keys=df_keys, axis=0).reset_index(level=1)
    # return hyper_scores_df, closeness_df, betweenness_df, rwr_df, fraction_betweenness_df
    return fraction_betweenness_df


def bridge_scores(adjacency_matrix, proteins_by_process_df, t_process):
    #proteins_by_process_df['protein_index'] = proteins_by_process_df['protein_index'].apply(literal_eval)
    t_process = t_process.replace('.', '-')
    protein_ppi_by_process = proteins_by_process_df.to_dict('list')
    totalppi_by_protein = list(np.sum(adjacency_matrix, axis=1))
    bridges = {}
    for process, ppis in protein_ppi_by_process.items():
        process = process.replace('.', '-')
        if process == t_process:
            continue
        else:
            proteins_bridge_by_process = [(ppi1*ppi2)/(totalppi) for ppi1, ppi2, totalppi in zip(
                protein_ppi_by_process[t_process], ppis, totalppi_by_protein)]
            bridges[process] = proteins_bridge_by_process
    t_process_bridge_df = pd.DataFrame.from_dict(
        bridges, orient='columns')
    return t_process_bridge_df


def scp(closeness_df, t_process):
    target_closeness = closeness_df[t_process].transpose().to_numpy()
    remaining_closeness = closeness_df.loc[:,
                                           closeness_df.columns != t_process].transpose().to_numpy()
    scp = target_closeness * remaining_closeness
    scp = pd.DataFrame(scp).transpose()
    #scp_df = pd.concat([scp_df, scp], axis=0).reset_index(drop=True)
    return scp


"""def slb(t_process, limit=5):
    graph = Graph.Read_GML("../python/data/graph_apid_huri")
    processes_df = pd.read_csv("../python/data/reactome_proteins_indexes_apid_huri.csv", sep=',', header=0)
    processes_df['protein_index'] = processes_df['protein_index'].apply(literal_eval)
    processes_df['proteins_ids'] = processes_df['proteins_ids'].apply(literal_eval)

    #graph = Graph.Read_GML("../data/graph_apid_huri")
    t_process = t_process.replace('.', '-')
    t_process_proteins = processes_df[processes_df['process']==t_process]['proteins_ids'].values[0]
    processes_df = processes_df[processes_df['process']!=t_process]
    remaining_process_proteins = list(set(
        [protein for protein_list in processes_df['proteins_ids'].values for protein in protein_list if protein in graph.vs['name']]))
    process_slb = {}
    for protein in tqdm(t_process_proteins):
        shortest_paths = graph.get_all_shortest_paths(protein, remaining_process_proteins)
        for index, process in processes_df.iterrows():
            process_id = process['process']
            if process_id not in process_slb:
                process_slb[process_id] = {}
            process_proteins = process['protein_index']
            process_shortest_paths = [shortest_path[1:-1] for shortest_path in shortest_paths if shortest_path[-1] in process_proteins and len(shortest_path)<=limit+1]
            for i in range(3,limit+1):
                shortest_paths_len = [shortest_path for shortest_path in process_shortest_paths if len(shortest_path)==i]
                shortest_path_i_count = Counter(chain.from_iterable(shortest_paths_len))
                shortest_path_i_count = {protein: n/i for protein,n in shortest_path_i_count.items()}
                process_slb[process_id] = dict(Counter(process_slb[process_id])+Counter(shortest_path_i_count))
            for protein in range(len(graph.vs['name'])):
                if protein not in process_slb[process_id]:
                    process_slb[process_id][protein] = 0
    process_slb_df =  pd.DataFrame.from_dict(process_slb, orient='columns')
    process_slb_df.sort_index(inplace=True)
    return process_slb_df"""


def slb(t_process, proteins_by_process=None, limit=5, reduced=False):
    if reduced:
        graph = Graph.Read_GML("../python/data/graph_apid_huri_80")
        processes_df = proteins_by_process
    else:
        graph = Graph.Read_GML("../data/graph_apid_huri")
        processes_df = pd.read_csv(
            "../data/reactome_proteins_indexes_apid_huri.csv", sep=',', header=0)
    try:
        processes_df['protein_index'] = processes_df['protein_index'].apply(
            literal_eval)
        processes_df['proteins_ids'] = processes_df['proteins_ids'].apply(
            literal_eval)
    except ValueError:
        processes_df['protein_index'] = processes_df['protein_index']
        processes_df['proteins_ids'] = processes_df['proteins_ids']
    t_process = t_process.replace('.', '-')
    print(t_process)
    if reduced:
        t_process_proteins = [protein for protein in processes_df[processes_df['process']
                                                                  == t_process]['proteins_ids'].values[0] if protein in graph.vs['name']]
    else:
        t_process_proteins = processes_df[processes_df['process']
                                          == t_process]['proteins_ids'].values[0]

    processes_df = processes_df[processes_df['process'] != t_process]
    protein_ids_processes_dict = {}
    for k, v in processes_df[['process', 'proteins_ids']].set_index('process').to_dict('dict')['proteins_ids'].items():
        for x in v:
            if x in graph.vs['name']:
                protein_ids_processes_dict.setdefault(x, []).append(k)

    protein_index_processes_dict = {}
    for k, v in protein_ids_processes_dict.items():
        protein_index_processes_dict[int(graph.vs.find(name=k)['id'])] = v

    process_slb = {}
    for process in processes_df['process']:
        process_slb[process] = {}

    shortest_path_dict = {}
    for protein in t_process_proteins:
        shortest_paths = graph.get_all_shortest_paths(
            protein, protein_ids_processes_dict.keys())
        for shortest_path in shortest_paths:
            if shortest_path[-1] not in shortest_path_dict:
                shortest_path_dict[shortest_path[-1]] = [shortest_path[1:-1]]
            else:
                shortest_path_dict[shortest_path[-1]
                                   ].append(shortest_path[1:-1])
    for end_protein, middle_proteins in shortest_path_dict.items():
        for i in range(1, limit):
            limit_middle_proteins = [
                middle_protein for middle_protein in middle_proteins if len(middle_protein) == i]
            shortest_path_i_count = Counter(
                chain.from_iterable(limit_middle_proteins))
            shortest_path_i_count = {
                protein: n/i for protein, n in shortest_path_i_count.items()}
            # print(protein_index_processes_dict)
            for process in protein_index_processes_dict[end_protein]:
                process_slb[process] = dict(
                    Counter(process_slb[process])+Counter(shortest_path_i_count))

    for process in processes_df['process']:
        for protein in range(len(graph.vs['name'])):
            if protein not in process_slb[process]:
                process_slb[process][protein] = 0

    process_slb_df = pd.DataFrame.from_dict(process_slb, orient='columns')
    process_slb_df.sort_index(inplace=True)
    return process_slb_df.to_numpy()


def slb2(t_process, proteins_by_process=None, limit=5, reduced=False):
    if reduced:
        graph = Graph.Read_GML("../python/data/graph_apid_huri_80")
        processes_df = proteins_by_process
        try:
            processes_df['protein_index'] = processes_df['protein_index'].apply(
                literal_eval)
            processes_df['proteins_ids'] = processes_df['proteins_ids'].apply(
                literal_eval)
        except ValueError:
            processes_df['protein_index'] = processes_df['protein_index']
            processes_df['proteins_ids'] = processes_df['proteins_ids']
        t_process_proteins = [protein for protein in processes_df[processes_df['process']
                                                                  == t_process]['proteins_ids'].values[0] if protein in graph.vs['name']]
    else:
        graph = Graph.Read_GML("../data/graph_apid_huri")
        processes_df = pd.read_csv(
            "../data/reactome_proteins_indexes_apid_huri.csv", sep=',', header=0)
        try:
            processes_df['protein_index'] = processes_df['protein_index'].apply(
                literal_eval)
            processes_df['proteins_ids'] = processes_df['proteins_ids'].apply(
                literal_eval)
        except ValueError:
            processes_df['protein_index'] = processes_df['protein_index']
            processes_df['proteins_ids'] = processes_df['proteins_ids']
        t_process_proteins = processes_df[processes_df['process']
                                          == t_process]['proteins_ids'].values[0]

    t_process = t_process.replace('.', '-')
    processes_df = processes_df[processes_df['process'] != t_process]
    protein_index_processes_dict = {}
    for k, v in processes_df[['process', 'proteins_ids']].set_index('process').to_dict('dict')['proteins_ids'].items():
        for x in v:
            if x in graph.vs['name']:
                protein_index_processes_dict.setdefault(
                    int(graph.vs.find(name=x)['id']), []).append(k)

    process_slb = {}
    for process in processes_df['process']:
        process_slb[process] = {}

    for protein in tqdm(t_process_proteins):

        shortest_paths = graph.get_all_shortest_paths(
            protein, protein_index_processes_dict.keys())
        for i in range(3, limit+2):
            length_shortest_paths = [
                short_path[1:] for short_path in shortest_paths if len(short_path) == i]
            shortest_path_dict = {}
            for length_shortest_path in length_shortest_paths:
                if length_shortest_path[-1] not in shortest_path_dict:
                    shortest_path_dict[length_shortest_path[-1]
                                       ] = length_shortest_path[:-1]
                else:
                    shortest_path_dict[length_shortest_path[-1]
                                       ].extend(length_shortest_path[:-1])
            for end_protein, middle_proteins in shortest_path_dict.items():
                shortest_path_i_count = Counter(middle_proteins)
                shortest_path_i_count = {
                    protein: n/i for protein, n in shortest_path_i_count.items()}
                for process in protein_index_processes_dict[end_protein]:
                    process_slb[process] = dict(
                        Counter(process_slb[process])+Counter(shortest_path_i_count))
    for process in processes_df['process']:
        for protein in range(len(graph.vs['name'])):
            if protein not in process_slb[process]:
                process_slb[process][protein] = 0

    process_slb_df = pd.DataFrame.from_dict(process_slb, orient='columns')
    process_slb_df.sort_index(inplace=True)
    return process_slb_df

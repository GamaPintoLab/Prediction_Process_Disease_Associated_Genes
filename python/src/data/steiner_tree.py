from igraph import Graph
import pandas as pd
from collections import Counter
import numpy as np


def steiner_tree(idlist, graph):
    starting_terminal = idlist[0]

    terminals_in_graph = [starting_terminal]
    tree_nodes = [starting_terminal]
    n = 0
    while any(terminal not in terminals_in_graph for terminal in idlist):
        # while n<10:
        subdist = graph.shortest_paths(source=tree_nodes, target=idlist)
        # print(subdist)
        for dist in subdist:
            for element in range(len(dist)):
                if dist[element] == 0:
                    dist[element] = 999

        min_dist = min([min(dist) for dist in subdist])
        new_node = -1
        while new_node == -1:
            print(min_dist)
            for term in range(len(subdist)):
                for node in range(len(subdist[term])):
                    # if subdist[term][node]==min_dist:
                    #print(term, node)
                    if idlist[node] not in terminals_in_graph and subdist[term][node] == min_dist:
                        new_node = idlist[node]
                        connection_node = term

                if new_node != -1:
                    break
            min_dist += 1
        if min_dist == 2:
            tree_nodes.append(new_node)
            terminals_in_graph.append(new_node)
        else:
            # print(connection_node)
            new_nodes = graph.get_all_shortest_paths(
                tree_nodes[connection_node], to=new_node)
            print(new_nodes)
            tree_nodes.extend(new_nodes[0][1:])
            terminals_in_graph.append(new_node)
            # print(terminals_in_graph)
        # print(len(terminals_in_graph))

        n += 1
        #shortest_paths = graph.get_all_shortest_paths(i, to=[terminal for terminal in idlist if terminal not in terminals_in_graph])
        # print(shortest_paths)
    print(len(tree_nodes))
    return


def sca(idlist, graph, adj):
    added_nodes = []
    increase = True
    original_nodes = idlist[:]
    while increase:
        subgraph = graph.induced_subgraph(idlist)
        components = Graph.clusters(subgraph, mode='strong')
        # print(list(components))
        comp_adj_matrix = adj[:, sorted(idlist)]
        lcc = max([len(component) for component in components])
        # print(lcc)
        ind_comp_adj_matrices = [comp_adj_matrix[:, comp]
                                 for comp in components]
        max_addition = []
        for ind_comp in ind_comp_adj_matrices:
            n_int = np.sum(ind_comp, axis=1)
            n_int[n_int >= 1] = ind_comp.shape[1]
            max_addition.append(n_int)
        max_addition = np.array(max_addition).transpose()
        max_addition_total = np.array(max_addition).sum(axis=1)
        if max(max_addition_total) > lcc:
            increase = True
            candidates = np.argwhere(max_addition_total == np.amax(
                max_addition_total)).flatten().tolist()
            if len(candidates) > 1:
                cand_dict = {}
                for cand in candidates:
                    cand_dict[cand] = len([id_ for id_ in idlist if id_ in np.argwhere(
                        adj[:, cand] == 1).flatten().tolist()])/len(np.argwhere(adj[:, cand] == 1).flatten().tolist())
                candidates = [key for key, value in cand_dict.items() if value == max(cand_dict.values())]

            idlist.append(candidates[0])
            added_nodes.append(candidates[0])

        else:
            increase = False
    subgraph = graph.induced_subgraph(idlist)
    final_components = list(Graph.clusters(subgraph, mode='strong'))
    final_component = max(final_components, key=len)
    main_component = [idlist[node] for node in final_component]
    conservative_module = list(set(original_nodes) & set(main_component))
    return pd.Series([main_component, conservative_module, added_nodes])


def sca_analytics(idlist, graph, adj, process):
    print()
    added_nodes = []
    increase = True
    original_nodes = idlist[:]
    df_dict = {'process':[], 'n_clusters': [], 'n_candidates':[], 'best_candidate':[], 'best_ratio':[]}
    while increase:
        subgraph = graph.induced_subgraph(idlist)
        components = Graph.clusters(subgraph, mode='strong')
        df_dict['process'].append(process)
        df_dict['n_clusters'].append(len(components))
        comp_adj_matrix = adj[:, sorted(idlist)]
        lcc = max([len(component) for component in components])
        ind_comp_adj_matrices = [comp_adj_matrix[:, comp]
                                 for comp in components]
        max_addition = []
        for ind_comp in ind_comp_adj_matrices:
            n_int = np.sum(ind_comp, axis=1)
            n_int[n_int >= 1] = ind_comp.shape[1]
            max_addition.append(n_int)
        max_addition = np.array(max_addition).transpose()
        max_addition_total = np.array(max_addition).sum(axis=1)
        if max(max_addition_total) > lcc:
            increase = True
            candidates = np.argwhere(max_addition_total == np.amax(
                max_addition_total)).flatten().tolist()
            df_dict['n_candidates'].append(len(candidates))
            cand_dict = {}
            for cand in candidates:
                cand_dict[cand] = len([id_ for id_ in idlist if id_ in np.argwhere(
                        adj[:, cand] == 1).flatten().tolist()])/len(np.argwhere(adj[:, cand] == 1).flatten().tolist())
            #df_dict['canditates_ratio'].append(cand_dict)
            candidates = [key for key, value in cand_dict.items() if value == max(cand_dict.values())]
            df_dict['best_candidate'].append(candidates[0])
            df_dict['best_ratio'].append(max(cand_dict.values()))
            idlist.append(candidates[0])
            added_nodes.append(candidates[0])

        else:
            increase = False
    df_dict['n_candidates'].append(0)
    df_dict['best_candidate'].append(0)
    df_dict['best_ratio'].append(0)
    subgraph = graph.induced_subgraph(idlist)
    final_components = list(Graph.clusters(subgraph, mode='strong'))
    final_component = max(final_components, key=len)
    main_component = [idlist[node] for node in final_component]
    conservative_module = list(set(original_nodes) & set(main_component))
    return pd.DataFrame.from_dict(df_dict)
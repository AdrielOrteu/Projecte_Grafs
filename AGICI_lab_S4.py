from Bio.SeqRecord import SeqRecord
from  Bio import SeqIO
import csv
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import math
import copy
from itertools import combinations
from random import sample
from networkx.algorithms.bipartite import color


# ------- IMPLEMENT HERE ANY AUXILIARY FUNCTIONS NEEDED ------- #


# --------------- END OF AUXILIARY FUNCTIONS ------------------ #

def feature_list(genome: SeqRecord, query: str) -> list:
    """
    Extract CDS features with specific feature description.
    - param genome : SeqRecord
        genome SeqRecord object to be analyzed.
    - param query : str
        feature descriptor.
    - return list
        list of tuples (locus_tag, protein_id) matching descriptor.
    """
    
    ret_list = []
    for feat in genome.features:
        if feat.type == 'CDS':
            if query in feat.qualifiers['product'][0]:
                if 'protein_id' in feat.qualifiers.keys():
                    ret_list.append((feat.qualifiers['locus_tag'], feat.qualifiers['protein_id']))
    return ret_list


def gene_qualifier(query: str, query_field: str,
                   target_field: str, genome: SeqRecord) -> tuple:
    """
    Obtain the specified qualifier identifier for a given gene qualifier
    - param: query: str
        gene name/locus_tag to map to corresponding locus_tag/name
    - param: query_field: str
        query type indicator (gene/locus_tag/protein_id/product)
    - param: target_field: str
        target type indicator (gen e/locus_tag/protein_id/product)
    - param: genome : SeqRecord
        genome SeqRecord object containing features
    - return: tuple
        int : feature index
        str : specified qualifier for gene (empty string if no match)
    """

    feat_num = 0
    ret_value = ''
    max_features = len(genome.features)
    while ret_value == '' and feat_num < max_features:
        if genome.features[feat_num].type == 'gene':
            if genome.features[feat_num].qualifiers[query_field][0] == query:
                ret_value = genome.features[feat_num].qualifiers[target_field][0]
            else:
                feat_num = feat_num + 1
        else:
            feat_num = feat_num + 1
    return feat_num, ret_value

'''
i = index
referencia = genome[index].final + 100
while genome[i].inical < referencia :
    llista.genome[i].locustag
    i += 1
'''
def operon(locus_tag: str='', max_intergenic_dist: int=100, genome: SeqRecord=None) -> list:
    '''
    Obtain the locus_tag identifier for a given gene name
    - param: locus_tag : str
        locus_tag of lead operon gene
    - param: max_intergenic_dist : int
        maximum distance between consecutive, same strand genes
    - param: genome : SeqRecord
        genome SeqRecord object for operon inference
    - return: list
        list of locus_tags conforming operon (including query)
    '''
    operon = []
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #

    index, locus_tag_value = gene_qualifier(query=locus_tag, query_field='locus_tag', target_field='locus_tag', genome=genome)
    if locus_tag_value == "":  # Comprovamos que se ha encontrado el gen/locus_tag que buscamos
        return operon
    gene_feat = genome.features[index]
    gene_strand = gene_feat.location.strand
    
    genes = [f for f in genome.features if f.type == "gene" and f.location.strand == gene_strand]
    genes.sort(key=lambda f: int(f.location.start), reverse=(gene_strand==-1) )
    
    start_index = genes.index(gene_feat)
    prev_pos = int(gene_feat.location.end) if gene_strand == 1 else int(gene_feat.location.start)
    for gene in genes[start_index+1:]:
        curr_pos = int(gene.location.start) if gene_strand == 1 else int(gene.location.end)
        if abs(curr_pos - prev_pos) <= max_intergenic_dist:
            operon.append(gene.qualifiers['locus_tag'][0])
            prev_pos = int(gene.location.end) if gene_strand == 1 else int(gene.location.start)
        else:
            break
    
    # ----------------- END OF FUNCTION --------------------- #
    return operon

def degree_distribution(G : nx.DiGraph, degree_type : str, \
        node_type : str, bestN : int) -> tuple:
        '''Compute the in-, out- or general degree distribution.
        Parameters
        ----------
        - param: G : Networkx graph
        Graph to analyze
        - param: degree_type : str
        Type of degree to compute ('in', 'out', 'general')
        - param: node_type : str
        Type of nodes on which we report degrees ('TG', 'TG', 'all')

        - return: list
        list with degree frequency distribution
        - return: dict
        dictionary with node names as keys and their degrees as values
        for the top N nodes in the degree distribution

        '''
        # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
        data_nodes = [nodes[0] for nodes in G.nodes(data=True)]
        if node_type == 'TF':
            data_nodes = [nodes[0] for nodes in [nodes for nodes in G.nodes(data=True)] if nodes[1]['type'] == 'TF']
        if node_type == 'TG':
            data_nodes = [nodes[0] for nodes in [nodes for nodes in G.nodes(data=True)] if nodes[1]['type'] == 'TG']
        if degree_type == 'in':
            degrees = [ x for x in G.in_degree() if x[0] in data_nodes]
        if degree_type == 'out':
            degrees = [ x for x in G.out_degree() if x[0] in data_nodes]
        if degree_type == 'general':
            degrees = [ x for x in G.degree() if x[0] in data_nodes]
        frequencias = Counter([x[1] for x in degrees])
        bestitems = sorted(degrees, key=lambda  x : x[1], reverse=True)[:bestN]
        return frequencias, bestitems
        # ----------------- END OF FUNCTION --------------------- #

def TF_RISet_parse(tf_riset_filename: str, tf_set_filename: str,
                   detect_operons: bool, max_intergenic_dist: int,
                   genome: SeqRecord, mode:bool = False) -> nx.DiGraph:
    """
    Parse TF-RISet file to obtain a TRN graph.
    The TFSet file will be used to extract information on the gene coding for
    each transcription factor [4)geneCodingForTF, 5)geneBnumberCodingForTF]
    The TF-RISet file will be used to extract genes regulated by each TF.
    We will create nodes using their locus_tag identifier, and save the gene
    name as a 'name' attribute.
    If selected, operons will be predicted for each of these genes to determine
    the entire set of genes regulated by the TF.

    - param: tf_riset_filename : str
        name/path of the TF_RISet file
    - param: tf_set_filename : str
        name/path of the TFSet file
    - param: detect_operons : bool
        whether we run operon detection
    - param: max_intergenic_dist : int
        maximum distance between consecutive, same strand genes
    - param: genome : SeqRecord
        genome SeqRecord object to extract information from
    - return: TF dictionary
    """
    # a = feature_list(genome=genome,query="a")
    if mode:
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    tf_dict = {}
    nodes = {}
    with open(tf_set_filename, newline='') as file: #creacio tf_dict
        tsv_reader = csv.reader(file, delimiter='\t')
        for index,row in enumerate(tsv_reader):
            if len(row) < 18:   
                continue
            if index != 35:
                tf_dict[row[1]] = row[4]


    with open (tf_riset_filename, newline='') as file:
        tsv_reader = csv.reader(file,delimiter='\t')
        for row in tsv_reader:
            if len(row) < 27:
                continue
            if index > 45 and row[3] in tf_dict:
                node1 = row[3]
                node2 = row[16]
                if node1 not in nodes:
                    nodes[node1] = True
                    G.add_node((tf_dict[node1]), color="red", name='node1', type='TF')
                TG_locus_tag = gene_qualifier(node2, 'gene', 'locus_tag', genome)[1]
                if detect_operons:
                    operon_genes = operon(TG_locus_tag,100,genome)
                if node2 not in nodes:
                    nodes[node2] = True
                    info_gene = feature_list(genome,node2)
                    G.add_node((TG_locus_tag), type='TG')
                G.add_edge(tf_dict[node1],TG_locus_tag)
                if detect_operons:
                    operones = (operon(locus_tag=TG_locus_tag, genome=genome))
                    if operones:
                        G.add_nodes_from(operones, type='TG')
                    for genes in operones:
                        G.add_edge(node1,genes)
                






    # ----------------- END OF FUNCTION --------------------- #

    return G

def largest_CC_graph(G : nx.Graph) -> nx.Graph:
 '''Generate the largest connected component graph.
 Parameters
 ----------
 - param: G : Networkx graph
 Graph to analyze
 - return: Networkx graph
 the graph corresponding to the largest connected component in G
 '''
 # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
 x = copy.deepcopy(G).to_undirected()
 #print(nx.number_connected_components(x))
 subG = sorted(nx.connected_components(x), key=len, reverse=True)[0]
 subG = x.subgraph(subG)
 #nx.draw(x)
 #print("Num edges de la CC més gran:", subG.number_of_edges())
 #print("Num nodes de la CC més gran:",subG.number_of_nodes())
 return subG
 # ----------------- END OF FUNCTION --------------------- #

def average_distance(G: nx.Graph, iterations: int) -> float:
    '''Estimate the average distance in the input graph.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - param: iterations : int
        Number of iterations to perform
    - return: float
        The estimated average distance in G
    '''

    avg_distance = 0
    node_list = list(G.nodes())
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    for i in range (0, iterations):
        random_nodes = sample(node_list,2)
        has_path = nx.has_path(G,random_nodes[0],random_nodes[1])
        while not has_path:
            random_nodes = sample(node_list,2)
            has_path = nx.has_path(G,random_nodes[0],random_nodes[1])
        avg_distance += len(nx.shortest_path(G, random_nodes[0], random_nodes[1]))-1
    # ----------------- END OF FUNCTION --------------------- #

    return avg_distance/iterations

def deletion_impact(G: nx.Graph, node_list: list, \
                    grouping_size: int, iterations: int) -> dict:
    '''Assess the impact of node deletions on the graph average distance.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - param: node_list : list
        List of nodes to delete from the network
    - param: grouping_size : list
        The size of the groupings among nodes in the list
    - param: iterations : int
        Number of iterations to perform for average distance
    - return: dict
        Dictionary with grouping node names tuples as keys and differential
        average distance as values.
    '''


    

    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    del_impact = {}
    per_quitar = combinations(node_list, grouping_size)
    ave_distance = nx.average_shortest_path_length(G)
    for group in per_quitar:
        nodes_finals = [x for x in G.nodes() if x not in per_quitar]
        sub_graf = G.subgraph(nodes_finals)
        del_impact[group] = ave_distance - average_distance(sub_graf, iterations)
    # ----------------- END OF FUNCTION --------------------- #
    return del_impact



if __name__ == "__main__":

    # ------- IMPLEMENT HERE THE MAIN FOR THIS SESSION ------- #
    import time
    start_time = time.time()
    genome = SeqIO.read('dataset/sequence.gb' , 'genbank')
    
    G_no = TF_RISet_parse('dataset/TF-RISet.tsv', 'dataset/TFSet.tsv', False, 100, genome, True)
    #G_so = TF_RISet_parse('dataset/TF-RISet.tsv', 'dataset/TFSet.tsv', True, 100, genome,True)
    #G1=nx.read_graphml('Ecoli_TRN.graphml')
    #G2=nx.read_graphml('Ecoli_operon_TRN.graphml')
    large_CC = largest_CC_graph(G_no)
    #print(nx.is_directed)
    #aproximació_500 = average_distance(large_CC,500)
    aproximació_1000 = average_distance(large_CC,1000)
    #aproximació_2000 = average_distance(large_CC,2000)
    #aproximació_10000 = average_distance(large_CC,10000)
    #aproximació_100000 = average_distance(large_CC,100000)

    valor_real = nx.average_shortest_path_length(large_CC)
    print(aproximació_1000, valor_real)

    ...

    print("--- %s seconds ---" % (time.time() - start_time))
    # ------------------- END OF MAIN ------------------------ #


from Bio.SeqRecord import SeqRecord
from  Bio import SeqIO
import csv
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import math
import copy

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
    operon = [locus_tag]
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
                    G.add_node((tf_dict[node1]), color="red")
                TG_locus_tag = gene_qualifier(node2, 'gene', 'locus_tag', genome)[1]
                if detect_operons:
                    operon_genes = operon(TG_locus_tag,100,genome)
                if node2 not in nodes:
                    nodes[node2] = True
                    info_gene = feature_list(genome,node2)
                    G.add_node((TG_locus_tag))
                G.add_edge(node1,node2)
                if detect_operons:
                    operones = (operon(locus_tag=TG_locus_tag, genome=genome))
                    if operones:
                        G.add_nodes_from(operones)
                    for genes in operones:
                        G.add_edge(node1,operon)
                






    # ----------------- END OF FUNCTION --------------------- #

    return G


if __name__ == "__main__": # temporarily disabled main to test specific parts of the code

    # ------- IMPLEMENT HERE THE MAIN FOR THIS SESSION ------- #

    import time
    start_time = time.time()

    # read in genome file as SeqRecord object
    genome = SeqIO.read('dataset/sequence.gb' , 'genbank')

    # parse TF_RISet file to obtain networks
    #miniTF-RISet.tsv para pruebas pequeñas
    G1 = TF_RISet_parse('dataset/TF-RISet.tsv', 'dataset/TFSet.tsv', False, 100, genome,True)
    Gpruebas = TF_RISet_parse('dataset/TF-RISet.tsv', 'dataset/TFSet.tsv', False, 100, genome,True) #Graf per mirar PowerLaw
    
    # report basic network stats
    '''
    Hay que añadir otras pruebas tipo grado medio camino mas corto medio y por delante
    '''
    #degress = [ x[1] for x in Gpruebas.degree] #Sagar el grado de cada nodo
    #frequencias = Counter(degress) # Sacar el Counter {valor  del grado: Cuantos veces sale el valor}
    #list_x = list(frequencias) # Sacar las claves del counter (valores del grado)
    #list_y = list(frequencias.values()) # Sacar los valores del counter (frequencias)
    #y_punts = np.array(copy.deepcopy(list_y)) # Valores para grafo de puntos normal x
    #x_punts = np.array(copy.deepcopy(list_x)) # Valores para grafo de puntos normal y
    #x_punts_plaw = np.array(list(map(math.log,copy.deepcopy(list_x)))) #Valores para power law x
    #y_punts_plaw = np.array(list(map(math.log, copy.deepcopy(list_y)))) #Valores para power law y
    #fig, ax = plt.subplots(2, 1)
    #ax[0].plot(x_punts, y_punts, 'o')
    #ax[1].plot(x_punts_plaw, y_punts_plaw, 'o')
    #plt.show()
    # export graph
    nx.write_graphml(G1, 'Ecoli_TRN.graphml')
    node_colors = [G1.nodes[n].get('color', 'blue') for n in G1.nodes()]
    nx.draw(G1, node_color=node_colors)
    plt.show()
    print(len(G1))
    print(G1.number_of_edges())

    print("--- %s seconds ---" % (time.time() - start_time))

    # ------------------- END OF MAIN ------------------------ #

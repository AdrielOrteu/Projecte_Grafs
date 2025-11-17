from Bio.SeqRecord import SeqRecord
from  Bio import SeqIO
import csv
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import math
import copy

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
        target type indicator (gene/locus_tag/protein_id/product)
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
        for index, row in tsv_reader:
            if len(row) < 27:
                continue
            if index > 45 and row[3] in tf_dict:
                node1 = row[3]
                node2 = row[16]
                if node1 not in nodes:
                    nodes[node1] = True
                    G.add_node((tf_dict[node1]))
                if node2 not in nodes:
                    nodes[node2] = True
                    info_gene = feature_list(genome,node2)
                    G.add_node(node2)
                G.add_edge(tf_dict[node1],node2)
    
    # ----------------- END OF FUNCTION --------------------- #
    
    return G
    

if __name__ == "__main__":

    # ------- IMPLEMENT HERE THE MAIN FOR THIS SESSION ------- #

    import time
    start_time = time.time()

    # read in genome file as SeqRecord object
    genome = SeqIO.read('dataset/sequence.gb' , 'genbank')

    # parse TF_RISet file to obtain networks
    G1 = TF_RISet_parse('dataset/TF-RISet.tsv', 'dataset/TFSet.tsv', False, 100, genome,False)
    Gpruebas = TF_RISet_parse('dataset/TF-RISet.tsv', 'dataset/TFSet.tsv', False, 100, genome,False) #Graf per mirar PowerLaw

    # report basic network stats
    '''
    Hay que añadir otras pruebas tipo grado medio camino mas corto medio y por delante
    '''
    degress = [ x[1] for x in Gpruebas.degree] #Sagar el grado de cada nodo
    frequencias = Counter(degress) # Sacar el Counter {valor  del grado: Cuantos veces sale el valor}
    list_x = list(frequencias) # Sacar las claves del counter (valores del grado)
    list_y = list(frequencias.values()) # Sacar los valores del counter (frequencias)
    y_punts = np.array(copy.deepcopy(list_y)) # Valores para grafo de puntos normal x
    x_punts = np.array(copy.deepcopy(list_x)) # Valores para grafo de puntos normal y
    x_punts_plaw = np.array(list(map(math.log,copy.deepcopy(list_x)))) #Valores para power law x
    y_punts_plaw = np.array(list(map(math.log, copy.deepcopy(list_y)))) #Valores para power law y
    fig, ax = plt.subplots(2, 1)
    ax[0].plot(x_punts, y_punts, 'o')
    ax[1].plot(x_punts_plaw, y_punts_plaw, 'o')
    plt.show()
    # export graph
    nx.write_graphml(G1, 'Ecoli_TRN.graphml')
    nx.draw_networkx(G1)

    print("--- %s seconds ---" % (time.time() - start_time))

    # ------------------- END OF MAIN ------------------------ #

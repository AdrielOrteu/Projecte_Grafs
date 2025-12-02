# -*- coding: utf-8 -*-

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import csv
import networkx as nx
import matplotlib.pyplot as plt


# ------- IMPLEMENT HERE ANY AUXILIARY FUNCTIONS NEEDED ------- #




# --------------- END OF AUXILIARY FUNCTIONS ------------------ #

def degree_distribution(G: nx.DiGraph, degree_type: str, \
                        node_type: str, bestN: int) -> tuple:
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

    deg_freq, bestN_degrees = [], {}

    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #





    # ----------------- END OF FUNCTION --------------------- #

    return deg_freq, bestN_degrees

def largest_CC_graph(G: nx.Graph) -> nx.Graph:
    '''Generate the largest connected component graph.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - return: Networkx graph
        The graph corresponding to the largest connected component in G
    '''

    H = None
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #





    # ----------------- END OF FUNCTION --------------------- #

    return H

if __name__ == "__main__":
    # ------- IMPLEMENT HERE THE MAIN FOR THIS SESSION ------- #
    import time
    start_time = time.time()

    G1=nx.read_graphml('Ecoli_TRN.graphml')
    G2=nx.read_graphml('Ecoli_operon_TRN.graphml')

    ...

    print("--- %s seconds ---" % (time.time() - start_time))
    # ------------------- END OF MAIN ------------------------ #


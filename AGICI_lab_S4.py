# -*- coding: utf-8 -*-

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import csv
import networkx as nx
import matplotlib.pyplot as plt


# ------- IMPLEMENT HERE ANY AUXILIARY FUNCTIONS NEEDED ------- #




# --------------- END OF AUXILIARY FUNCTIONS ------------------ #

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

    avg_distance = None
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #





    # ----------------- END OF FUNCTION --------------------- #

    return avg_distance

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

    del_impact = {}
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #






    # ----------------- END OF FUNCTION --------------------- #
    return del_impact

if __name__ == "__main__":

    # ------- IMPLEMENT HERE THE MAIN FOR THIS SESSION ------- #
    import time
    start_time = time.time()

    G1=nx.read_graphml('Ecoli_TRN.graphml')
    G2=nx.read_graphml('Ecoli_operon_TRN.graphml')

    ...

    print("--- %s seconds ---" % (time.time() - start_time))
    # ------------------- END OF MAIN ------------------------ #


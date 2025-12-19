import sys
import networkx as nx
import numpy as np
from itertools import combinations
from typing import List, Tuple
from io import StringIO

def find_motifs_FFL(G: nx.DiGraph) -> List[Tuple[str, str, str, str]]:
    """
    Detect motifs in the directed graph given as input.

    Parameters
    ----------
    - param: G : Networkx digraph
        Graph to analyze
    - return: list
        List of detected motifs, each represented as a tuple (X, Y, Z).
    """

    motifs = []
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    for nodes in G.nodes():
        for neighbor in G.successors(nodes):
            if neighbor != nodes:
                for neighbor_3 in G.successors(neighbor):
                    if neighbor_3 != neighbor and neighbor_3 != nodes:
                        if neighbor_3 in G.successors(nodes):
                            motifs.append((nodes,neighbor,neighbor_3))
    
    # ----------------- END OF FUNCTION --------------------- #
    return motifs

if __name__ == "__main__":

    # Read all stdin as the GraphML content
    # E.g1.: python3 AGICI_lab_S5.py < toy_1.graphml
    graphml_content = sys.stdin.read()
    # Load GraphML content a networkx graph
    G = nx.read_graphml(StringIO(graphml_content))
    #G = nx.read_graphml("mini_Ecoli_TRN.graphml")


    # Print number of motifs
    motifs = find_motifs(G)


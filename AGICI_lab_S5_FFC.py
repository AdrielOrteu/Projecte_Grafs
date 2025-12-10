import sys
import networkx as nx
import numpy as np
from itertools import combinations
from typing import List, Tuple
from io import StringIO

def find_motifs(G: nx.DiGraph) -> List[Tuple[str, str, str, str]]:
    """
    Detect motifs in the directed graph given as input.

    Parameters
    ----------
    - param: G : Networkx digraph
        Graph to analyze
    - return: list
        List of detected motifs, each represented as a tuple (X, Y, Z).
    """

def find_motifs(G: nx.DiGraph) -> List[Tuple[str, str, str, str]]:
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
<<<<<<< HEAD:AGICI_lab_BF.py
    G_nodes = {nom for (nom,val) in G.out_degree() if val >= 2}

    for n1, n2 in combinations(G_nodes, 2):
        intersect = set({x for x in set(G.successors(n1)) if x != n1 and x != n2}) & set({x for x in set(G.successors(n2)) if x != n1 and x != n2})
        for n3, n4 in combinations(intersect, 2):
            motifs.append((n1,n2,n3,n4))
=======
    for nodes in G.nodes():
        for neighbor in G.successors(nodes):
            if neighbor != nodes:
                for neighbor_3 in G.successors(neighbor):
                    if neighbor_3 != neighbor and neighbor_3 != nodes:
                        if neighbor_3 in G.successors(nodes):
                            motifs.append((nodes,neighbor,neighbor_3))
    
>>>>>>> b6595fa8c41aaa59690fb441b7788d81b23ed653:AGICI_lab_S5_FFC.py
    # ----------------- END OF FUNCTION --------------------- #
    print(len(motifs))
    return motifs

if __name__ == "__main__":

    # Read all stdin as the GraphML content
    # E.g.: python3 AGICI_lab_S5.py < toy_1.graphml
    graphml_content = sys.stdin.read()
    # Load GraphML content a networkx graph
    G = nx.read_graphml(StringIO(graphml_content))
<<<<<<< HEAD:AGICI_lab_BF.py
=======
    #G = nx.read_graphml("mini_Ecoli_TRN.graphml")


>>>>>>> b6595fa8c41aaa59690fb441b7788d81b23ed653:AGICI_lab_S5_FFC.py
    # Print number of motifs
    motifs = find_motifs(G)


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

param: G : Networkx digraph
      Graph to analyze
return: list
    List of detected motifs, each represented as a tuple (X, Y, Z)."""
    
    motifs = []
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    
    G_nodes = [nom for (nom, val) in G.out_degree() if val >= 2]
    
    for n1, n2 in combinations(G_nodes, 2):
        intersect = set(G.successors(n1)) & set(G.successors(n2))
        for n3, n4 in combinations(intersect, 2):
            motifs.append((n1, n2, n3, n4))
    # ----------------- END OF FUNCTION --------------------- #
    print(len(motifs))
    return motifs


if __name__ == "__main__":
    g1 = nx.DiGraph()
    g1.add_edge(1, 2)
    g1.add_edge(1, 3)
    g1.add_edge(1, 4)
    g1.add_edge(4, 2)
    g1.add_edge(4, 3)
    g1.add_edge(2, 3)
    g1.add_edge(2, 4)
    
    g2 = None
    
    # Read all stdin as the GraphML content
    # E.g1.: python3 AGICI_lab_S5.py < toy_1.graphml
    #graphml_content = sys.stdin.read()
    # Load GraphML content a networkx graph
    #G = nx.read_graphml(StringIO(graphml_content))
    # Print number of motifs
    motifs = find_motifs(g1)
    a = {1,2,3}
    b = {-3, -2, -1}
    print(tuple(a | b))

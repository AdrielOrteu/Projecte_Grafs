import sys
import networkx as nx
import numpy as np

from typing import List, Tuple
from io import StringIO

def find_motifs(G: nx.DiGraph) -> List[Tuple[str, str, str]]:
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
    
    
    # ----------------- END OF FUNCTION --------------------- #
    print(len(motifs))
    return motifs

if __name__ == "__main__":

    # Read all stdin as the GraphML content
    # E.g1.: python3 AGICI_lab_S5.py < toy_1.graphml
    graphml_content = sys.stdin.read()

    # Load GraphML content a networkx graph
    G = nx.read_graphml(StringIO(graphml_content))

    # Print number of motifs
    motifs = find_motifs(G)

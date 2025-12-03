import networkx as nx
import matplotlib.pyplot as plt

# Create a graph
G = nx.Graph()

# Define your node groups
nodes1 = [1, 2, 3]
nodes2 = [4, 5]

# Add nodes with colors
G.add_nodes_from(nodes1, color='red')
G.add_nodes_from(nodes2, color='green')

# Optional: add some edges
G.add_edges_from([(1,4), (2,5), (3,4)])

# Get the colors from node attributes
node_colors = [G.nodes[n].get('color', 'blue') for n in G.nodes()]

# Draw the graph
nx.draw(G, with_labels=True, node_color=node_colors)
plt.show()

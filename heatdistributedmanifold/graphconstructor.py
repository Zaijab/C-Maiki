import networkx as nx
import numpy as np
import string

dt = [('len', float)]
A = np.load("/home/zjabbar/code/C-Maiki/heatdistributedmanifold/unifracdist.npy")
A = A.view(dt)

G = nx.from_numpy_matrix(A)
G = nx.relabel_nodes(G, dict(zip(range(len(G.nodes())),string.ascii_uppercase)))    

G = nx.drawing.nx_agraph.to_agraph(G)

G.node_attr.update(color="red", style="filled")
G.edge_attr.update(color="blue", width="2.0")

G.draw('/tmp/out.png', format='png', prog='neato')

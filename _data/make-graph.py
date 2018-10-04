#!/usr/bin/env python

"""
    make-graph.py
    
    !! Graphs generated in this way aren't working ATM
"""

# from __future__ import print_function

# import sys
# import numpy as np
# import networkx as nx

# num_nodes = 100
# p = 0.10

# g = nx.erdos_renyi_graph(num_nodes, p=p)

# for node in nx.isolates(g):
#     g.remove_node(node)

# adj = nx.adjacency_matrix(g).tocoo()

# edgelist = np.column_stack([
#     adj.row,
#     adj.col,
# ])
# num_edges = edgelist.shape[0]

# with open('_data/edges.txt', 'w') as f:
#     f.write('%d %d\n' % (num_nodes, num_edges))
#     for src, dst in edgelist:
#         f.write('%d %d\n' % (src, dst))

# x = np.random.uniform(-100, 100, size=num_nodes)
# x[np.radom.choice(x.shape[0] * 0.5)] = 0
# np.savetxt('_data/nodes.txt', x, '%.3f')

# print('make-graph.py: wrote _data/{edges,nodes}.txt', file=sys.stderr)
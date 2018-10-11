#!/usr/bin/env python

"""
    make-graph.py
"""

from __future__ import print_function

import sys
import numpy as np
import networkx as nx

num_nodes = 100
p = 0.10

g = nx.erdos_renyi_graph(num_nodes, p=p, directed=True)

for node in nx.isolates(g):
    g.remove_node(node)

adj = nx.adjacency_matrix(g).tocoo()
adj.eliminate_zeros()

edgelist = np.column_stack([
    adj.row,
    adj.col,
])
num_edges = edgelist.shape[0]

with open('_data/edges.txt', 'w') as f:
    f.write('%d %d\n' % (num_nodes, num_edges))
    for src, dst in edgelist:
        f.write('%d %d\n' % (src, dst))

x = np.random.uniform(-10, 10, size=num_nodes)
np.savetxt('_data/nodes.txt', x, '%.3f')

print('make-graph.py: wrote _data/{edges,nodes}.txt', file=sys.stderr)
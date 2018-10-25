#!/usr/bin/env python

"""
    make-graph.py
"""

from __future__ import print_function

import sys
import numpy as np
import pandas as pd
import networkx as nx

num_nodes = 20
p = 0.5

g = nx.erdos_renyi_graph(num_nodes, p=p)

for node in list(nx.isolates(g)):
    g.remove_node(node)

adj = nx.adjacency_matrix(g).tocoo()

edgelist = np.column_stack([
    adj.row,
    adj.col,
])

num_nodes = adj.shape[0]

sel = edgelist[:,0] > edgelist[:,1]
edgelist[sel] = np.column_stack([edgelist[sel,1], edgelist[sel,0]])

edgelist = pd.DataFrame(edgelist).sort_values([0, 1]).drop_duplicates().values
num_edges = edgelist.shape[0]
num_nodes = np.unique(np.hstack(edgelist)).shape[0]

assert edgelist.max() == num_nodes - 1

with open('_data/edges.txt', 'w') as f:
    f.write('%d %d\n' % (num_nodes, num_edges))
    for src, dst in edgelist:
        f.write('%d %d\n' % (src, dst))

x = np.random.uniform(-100, 100, size=num_nodes)#.astype(int)
np.savetxt('_data/nodes.txt', x, '%.3f')

print('make-graph.py: wrote _data/{edges,nodes}.txt', file=sys.stderr)

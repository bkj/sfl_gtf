#!/usr/bin/env python

"""
    make-graph.py
"""

from __future__ import print_function

import sys
import numpy as np
import networkx as nx
import pandas as pd

num_nodes = 300
p = 0.01

g = nx.erdos_renyi_graph(num_nodes, p=p)

adj = nx.adjacency_matrix(g)
adj = ((adj + adj.T) > 0).astype(int)

keep = np.asarray(adj.sum(axis=-1)).squeeze() != 0
adj = adj[keep][:,keep]

adj = adj.tocoo()

edgelist = np.column_stack([
    adj.row,
    adj.col,
])

edgelist = edgelist[edgelist[:,0] < edgelist[:,1]]

num_edges = edgelist.shape[0]
num_nodes = adj.shape[0]

with open('_data/e', 'w') as f:
    f.write('%d %d\n' % (num_nodes, num_edges))
    for src, dst in edgelist:
        f.write('%d %d\n' % (src, dst))

x = np.random.uniform(-100, 100, size=num_nodes)

np.savetxt('_data/n', x, '%.3f')

print('make-graph.py: wrote _data/{e,n}', file=sys.stderr)

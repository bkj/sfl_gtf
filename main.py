import pandas as pd
import numpy as np
import networkx as nx

from heapq import heappush, heappop

np.set_printoptions(linewidth=130)

def prune(g):
    for edge in list(g.edges):
        if g.edges[edge]['weight'] == 0:
            g.remove_edge(*edge)
    
    return g

alpha = 3

def create_graph():
    edges = pd.read_csv('./_data/e', sep=' ', header=None, skiprows=1).values
    nodes = pd.read_csv('./_data/n', header=None).values.squeeze()
    assert len(set(np.hstack(edges))) == nodes.shape[0]
    
    num_nodes = nodes.shape[0]
    
    edges = np.vstack([
        edges,
        np.column_stack([edges[:,1], edges[:,0]])
    ])
    
    g = nx.from_edgelist(edges, create_using=nx.DiGraph())
    
    for edge in g.edges:
        g.edges[edge]['weight'] = alpha
        
    g.add_node('src')
    g.add_node('trg')
    
    node_mean = nodes.mean()
    
    for node_idx in range(num_nodes):
        if nodes[node_idx] > node_mean:
            g.add_edge('src', node_idx)
            g.edges[('src', node_idx)]['weight'] = nodes[node_idx] - node_mean
        elif nodes[node_idx] < node_mean:
            g.add_edge(node_idx, 'trg')
            g.edges[(node_idx, 'trg')]['weight'] = node_mean - nodes[node_idx]
        else:
            raise Exception
    
    return g, node_mean


def remove_edge(g, edge):
    if edge in g.edges:
        g.remove_edge(*edge)


def ensure_edge(g, edge):
    if edge not in g.edges:
        g.add_edge(*edge)
        g.edges[edge]['weight'] = 0


def bfs(g, parents):
    seen = set([])
    q = ['src']
    while len(q) > 0:
        u, q = q[0], q[1:]
        for v in sorted(g.neighbors(u)):
            if v not in seen and g.edges[(u, v)]['weight'] > 0:
                q.append(v)
                seen.add(v)
                parents[v] = u
    
    return 'trg' in seen


def max_flow(g):
    rg = g.copy()
    
    mf = 0
    counter = 0
    while True:
        parents = {}
        has_path = bfs(rg, parents)
        if not has_path:
            break
        
        # Find bottleneck
        path_flow = np.inf
        edge_path = []
        v = 'trg'
        while True:
            u = parents[v]
            path_flow = min(rg.edges[(u, v)]['weight'], path_flow)
            edge_path.append((u, v))
            if u == 'src':
                break
            else:
                v = u
        
        for (u, v) in edge_path:
            ensure_edge(rg, (u, v))
            ensure_edge(rg, (v, u))
            rg.edges[(u, v)]['weight'] -= path_flow
            rg.edges[(v, u)]['weight'] += path_flow
        
        counter += 1
        mf += path_flow
    
    rg = prune(rg)
    return rg, mf


g, node_mean = create_graph()

from collections import defaultdict

nlab = 1
averages      = defaultdict(lambda : 0)
nextlabel     = defaultdict(lambda : 0)
nums          = defaultdict(lambda : 0)
values        = defaultdict(lambda : 0)
inactivelabel = defaultdict(lambda : 0)
values        = {0 : node_mean}
label         = np.zeros(len(g.nodes), dtype=int)
alive         = dict(zip(g.nodes, [True] * len(g.nodes)))

it = 0
flagstop = 1
while flagstop:
    g, mf = max_flow(g)
    visited = nx.descendants(g, 'src')
    
    for n in range(nlab):
        averages[n] = 0
        nums[n] = 0
        nextlabel[n] = 0
    
    oldnlab = nlab
    for node in g.nodes:
        if node != 'src' and node != 'trg':
            if alive[node]:
                if node in visited:
                    curr_lab = label[node]
                    l = nextlabel[curr_lab]
                    if l == 0:
                        l = nlab
                        nextlabel[curr_lab] = nlab
                        nlab        += 1
                        averages[l]      = 0
                        nums[l]          = 0
                        nextlabel[l]     = 0
                        inactivelabel[l] = 0
                        values[l]        = values[curr_lab]
                    
                    label[node] = l
                    if ('src', node) in g.edges:
                        averages[l] += g.edges[('src', node)]['weight']
                    nums[l] += 1
                else:
                    l = label[node]
                    if (node, 'trg') in g.edges:
                        averages[l] -= g.edges[(node, 'trg')]['weight']
                    
                    for n in list(g.neighbors(node)):
                        if n in visited:
                            remove_edge(g, (node, n))
                    
                    nums[l] += 1

    for l in nums.keys():
        if l < oldnlab:
            if not inactivelabel[l]:
                if nextlabel[l] == 0:
                    averages[l] = 0
                    inactivelabel[l] = 1
                elif nums[l] == 0:
                    inactivelabel[l] = 1
                    inactivelabel[nextlabel[l]] = 1
                    averages[nextlabel[l]] = 0
                else:
                    averages[l] = averages[l] / nums[l]
                    values[l] += averages[l]
            else:
                averages[l] = 0
        else:
            averages[l] = averages[l] / nums[l]
            values[l] += averages[l]

    flagstop = 0
    for node in g.nodes:
        if node != 'src' and node != 'trg':
            if alive[node]:
                l = label[node]
                if inactivelabel[l]:
                    if node in visited:
                        remove_edge(g, ('src', node))
                    else:
                        remove_edge(g, (node, 'trg'))
                
                    alive[node] = 0
                    inactivelabel[l] = 1
                else:
                    flagstop = 1
                    if node in visited:
                        ensure_edge(g, ('src', node))
                        g.edges[('src', node)]['weight'] -= averages[l]
                        if g.edges[('src', node)]['weight'] < 0:
                            tmp = g.edges[('src', node)]['weight']
                            ensure_edge(g, (node, 'trg'))
                            g.edges[('src', node)]['weight'] = g.edges[(node, 'trg')]['weight']
                            g.edges[(node, 'trg')]['weight'] = -tmp
                    else:
                        ensure_edge(g, (node, 'trg'))
                        g.edges[(node, 'trg')]['weight'] += averages[l]
                        if g.edges[(node, 'trg')]['weight'] < 0:
                            tmp = g.edges[(node, 'trg')]['weight']
                            ensure_edge(g, ('src', node))
                            g.edges[(node, 'trg')]['weight'] = g.edges[('src', node)]['weight']
                            g.edges[('src', node)]['weight'] = -tmp
                            
    prune(g)
    
    it += 1

for node in g.nodes:
    if node != 'src' and node != 'trg':
        val = values[label[node]]
        val = max(val - alpha, 0) + min(val + alpha, 0)
        print(val)


#!/usr/bin/env python

"""
    sfl.py
"""

from __future__ import print_function

import sys
import argparse
import numpy as np
import pandas as pd

from supervx import SuperVX, SuperGraph

# --
# Args

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--node-path', type=str)
    parser.add_argument('--edge-path', type=str)
    parser.add_argument('--reg-sparsity', type=float, default=3)  # 3
    parser.add_argument('--reg-edge', type=float, default=3)      # 3
    parser.add_argument('--validate', type=str)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    
    edges = pd.read_csv(args.edge_path, header=None, sep=' ', skiprows=1).values
    feats = np.array([float(xx) for xx in open(args.node_path).read().splitlines()])
    
    supergraph = SuperGraph(
        edges=edges,
        feats=feats, 
        partition_mode=None,
        num_partitions=1,
    )
    
    supervx = SuperVX(
        supernodes=supergraph.supernodes,
        superedges=supergraph.superedges,
        reg_sparsity=args.reg_sparsity,
        reg_edge=args.reg_edge,
    )
    
    if not args.validate:
        _ = supervx.solve(UseADMM=False, Verbose=False, Rho=1e-10, EpsAbs=1e-12, EpsRel=1e-12)
        fitted = supergraph.unpack(supervx.values)
        fitted[np.abs(fitted) < 0.01] = 0
        print('\n'.join(fitted.astype(str)))
        
        obj = supervx.gvx.GetTotalProblemValue()
    else:
        fitted = pd.read_csv(args.validate, header=None).values.squeeze()
        sfeats = supergraph.supernodes.items()[0][1]
        sedges = np.vstack(supergraph.superedges[(1, 1)])
        obj = (
            0.5 * ((fitted - sfeats) ** 2).sum() + 
            args.reg_sparsity * np.abs(fitted).sum() + 
            args.reg_edge * np.abs(fitted[sedges[:,0]] - fitted[sedges[:,1]]).sum()
        )
    
    print("objective=%f" % obj, file=sys.stderr)
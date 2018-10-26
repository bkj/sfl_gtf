#!/usr/bin/env python

"""
    validate.py
"""

from __future__ import print_function

import sys
import numpy as np

a = np.array([float(xx) for xx in open(sys.argv[1]).read().splitlines()])
b = np.array([float(xx) for xx in open(sys.argv[2]).read().splitlines()])

a = a.round(3)
b = b.round(3)

z = np.column_stack([a, b])
num_nnz = (z.sum(axis=-1) != 0).sum()
print('nonzero entries: %d' % num_nnz)
# print(z[z.sum(axis=-1) != 0])

if np.allclose(a, b):
    print("PASS")
else:
    print("FAIL")
    num_err = 0
    for aa, bb in zip(a, b):
        if not np.allclose([aa], [bb]):
            num_err += 1
            if num_err <= 10:
                print('%f != %f' % (aa, bb), file=sys.stderr)

print('num_err=%d' % num_err)

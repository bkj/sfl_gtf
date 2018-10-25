#!/usr/bin/env python

"""
    validate.py
"""

from __future__ import print_function

import sys
import numpy as np

a = np.array([float(xx) for xx in open(sys.argv[1]).read().splitlines()])
b = np.array([float(xx) for xx in open(sys.argv[2]).read().splitlines()])

# a = a.round(2)
# b = b.round(2)

z = np.column_stack([a, b])
print('nonzero entries: %d' % (z.sum(axis=-1) != 0).sum())
# print(z[z.sum(axis=-1) != 0])

if np.allclose(a, b):
    print("PASS")
else:
    print("FAIL")
    max_diff = np.max(np.abs(a - b))
    argmax_diff = np.argmax(np.abs(a - b))
    print(max_diff)
    print(a[argmax_diff], b[argmax_diff])
    # num_err = 0
    # for aa, bb in zip(a, b):
    #     if not np.allclose([aa], [bb]):
    #         num_err += 1
    #         print('%f != %f' % (aa, bb), file=sys.stderr)

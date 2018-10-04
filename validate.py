#!/usr/bin/env python

"""
    validate.py
"""

from __future__ import print_function

import sys
import numpy as np

a = [float(xx) for xx in open(sys.argv[1]).read().splitlines()]
b = [float(xx) for xx in open(sys.argv[2]).read().splitlines()]

if np.allclose(a, b):
    print("PASS")
else:
    print("FAIL")
    num_err = 0
    for aa, bb in zip(a, b):
        if not np.allclose([aa], [bb]):
            num_err += 1
            print('%f != %f' % (aa, bb), file=sys.stderr)

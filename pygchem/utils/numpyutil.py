# -*- coding: utf-8 -*-

# module numpy_addon
# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013 Benoît Bovy
# see license.txt for more details
# 
#
# Last modification: 06/2013

"""
Miscellaneous routine(s) based on the numpy package
"""

import numpy as np


def broadcast_1d_array(arr, ndim, axis=1):
    """
    Broadcast 1-d array `arr` to `ndim` dimensions on the first axis
    (`axis`=0) or on the last axis (`axis`=1).
     
    Useful for 'outer' calculations involving 1-d arrays that are related to
    different axes on a multidimensional grid.
    """
    ext_arr = arr
    for i in xrange(ndim - 1):
        ext_arr = np.expand_dims(ext_arr, axis=axis)
    return ext_arr

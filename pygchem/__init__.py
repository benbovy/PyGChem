# -*- coding: utf-8 -*-
"""
PyGChem
=======

This is PyGChem, a python interface to the
GEOS-Chem Global Chemistry Model.

It provides a high-level, Object-Oriented interface to the model
for input/output manipulation, visualizing, pre and post-processing.

pygchem is written by Benoît Bovy (bbovy at ulg.ac.be).
Code written by Gerrit Kuhlmann is also included in this package.
:copyright: Copyright 2013 by Benoît Bovy.
:license: GPL, see LICENSE for details.

"""

__all__ = ["globchem", "diagnostics", "io", "utils", ]
__version__ = (0, 0, 1)

# Check python version
import sys
if sys.version_info[:2] < (2, 6):
    raise ImportError("Python version 2.6 or later is required for PyGChem"
                      " (%d.%d detected)." % sys.version_info[:2])

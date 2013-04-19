# -*- coding: utf-8 -*-

# module pygchem.utils.error_handling
# pygchem: Python interface for GEOS-Chem Chemistry Transport Model
#
# Copyright (C) 2012 Benoit Bovy
# see license.txt for more details
#
# Last modification: 01/2013

"""
Module for handling errors and warnings
"""

import sys
import warnings

def warnings_stdout(message, category, filename, lineno, file=None):
    sys.stdout.write('%s: %s\n' % (category.__name__, message))
    return

warnings.showwarning = warnings_stdout
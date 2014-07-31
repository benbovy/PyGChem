# -*- coding: utf-8 -*-

# module config
# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2012/2013 Benoit Bovy.
# see license.txt for more details
# 
#
# Last modification: 04/2013

"""
Package configuration. "Global constants" are defined here.

"""

import inspect
import os


# Paths
_this_file = inspect.getframeinfo(inspect.currentframe()).filename
PACKAGE_ROOT_PATH = os.path.dirname(os.path.abspath(_this_file))
PACKAGE_DATA_PATH = os.path.join(PACKAGE_ROOT_PATH, "data")
del _this_file

# Default GEOS-Chem version
GCHEM_VERSION = "v9-02"

# physical or chemical constants
C_MOLECULAR_WEIGHT = 12e-3    # molecular weight of C atoms (kg/mole)

# other
WARN_STDOUT = False
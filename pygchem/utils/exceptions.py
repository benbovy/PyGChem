# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013-2014 Benoît Bovy
# see license.txt for more details
#

"""
Exceptions specific to the PyGChem package.

"""

import warnings
import sys

from pygchem import config


class PyGChemError(Exception):
    """Base class for errors in the PyGChem package."""
    pass


class NotYetImplementedError(PyGChemError):
    """
    Raised by missing functionality.

    Different meaning to NotImplementedError, which is for abstract methods.

    """
    pass


class NotPermittedError(PyGChemError):
    """Raised when an operation is not permitted."""
    pass


class SelectionMismatchError(PyGChemError):
    """
    Raised when a selection operation has failed to find the correct number
    of results.

    """
    pass


def warnings_stdout(message, category, filename, lineno, file=None):
    sys.stdout.write('%s: %s\n' % (category.__name__, message))
    return

if config.WARN_STDOUT:
    warnings.showwarning = warnings_stdout


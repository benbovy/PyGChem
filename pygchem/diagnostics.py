# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013-2014 Benoît Bovy
# see license.txt for more details
#

"""
GEOS-Chem diagnostics (set/get metadata to/from files 'tracerinfo.dat'
and 'diaginfo.dat').

"""

import os

from pygchem import config
from pygchem.io import diaginfo
from pygchem.utils.data_structures import record_cls, RecordList
from pygchem.utils.exceptions import NotYetImplementedError


CTMCategory = record_cls(
    "CTMCategory",
    "A category for CTM diagnostics, as stored in diaginfo.dat",
    (('offset', int, 0, True, 'Offset (constant to add to tracer numbers '
                              'in order to distinguish for the given diagnostic'
                              ' category, as stored in file tracerinfo.dat)'),
     ('name', str, None, True, 'Name of the category'),
     ('description', str, '', True, 'Descriptive comment')),
    required_fields=('offset', 'name')
)


CTMDiagnostic = record_cls(
    "CTMDiagnostic",
    "A CTM diagnostic, as stored in tracerinfo.dat",
    (('number', int, None, True, 'Diagnostic number'),
     ('name', str, None, True, 'Name of the CTM diagnostic (tracer)'),
     ('full_name', str, '', True, 'Long name of the diagnostic'),
     ('unit', str, 'unitless', True, 'Diagnostic unit'),
     ('scale', float, 1., True, 'Standard scale factor to convert values to '
                                'the diagnostic\'s unit'),
     ('chemical', bool, True, True, 'True if the diagnostic is a chemical '
                                    'tracer'),
     ('molecular_weight', float, 0., True, 'Molecular weight for chemical '
                                           'tracer (kg/mole)'),
     ('hydrocarbon', bool, False, True, 'True if the diagnostic is an '
                                        'hydrocarbon tracer'),
     ('carbon_weight', int, 0, True, 'Number of moles C/moles tracer for '
                                     'hydrocarbon tracers')
     ),
    required_fields=('number', 'name')
)


class CTMDiagnosticInfo(object):
    """
    Define all the diagnostics (and its metadata) that one can archive
    with GEOS-Chem.

    An instance of this class is commonly related to a couple
    of diaginfo.dat and tracerinfo.dat files. The class provides methods to
    read/write information from/to these files.

    Parameters
    ----------
    diaginfo_file : string or None
        path to the 'diaginfo.dat' file.
        If None, no file read (no category added)
        If empty string, the instance is created using a default
        diaginfo.dat file - located in the data folder of the module
        installation path)
    tracerinfo_file : string or None
        path to the 'tracerinfo.dat' file (or None or empty string).

    """
    def __init__(self, diaginfo_file='', tracerinfo_file=''):

        self.categories = RecordList([], key_attr='name',
                                     ref_classes=CTMCategory)
        self.diagnostics = RecordList([], key_attr='number',
                                      ref_classes=CTMDiagnostic)

        if diaginfo_file is not None:
            if not diaginfo_file:
                diaginfo_file = os.path.join(config.PACKAGE_DATA_PATH,
                                             "diaginfo.dat")
            self.load_diaginfo(diaginfo_file)

        if tracerinfo_file is not None:
            if not tracerinfo_file:
                tracerinfo_file = os.path.join(config.PACKAGE_DATA_PATH,
                                               "tracerinfo.dat")
            self.load_tracerinfo(tracerinfo_file)

    def load_diaginfo(self, filename, clear=True):
        """
        load diagnostic categories metadata from the 'diaginfo.dat' file given
        by `filename`. If `clear` is True, all existing categories are removed.
        """
        data = diaginfo.read_diaginfo(filename)
        if clear:
            del self.categories[:]
        self.categories.extend(CTMCategory(**d) for d in data)

    def load_tracerinfo(self, filename, clear=True):
        """
        Read diagnostics from the 'tracerinfo.dat'-file with given `filename`.
        If `clear` is True, all existing diagnostics are removed.
        """
        data = diaginfo.read_tracerinfo(filename)
        if clear:
            del self.diagnostics[:]
        for d in data:
            if d['carbon_weight'] != 1:
                d['hydrocarbon'] = True
                d['molecular_weight'] = config.C_MOLECULAR_WEIGHT
            else:
                d['hydrocarbon'] = False
            d['chemical'] = bool(d['molecular_weight'])

        self.diagnostics.extend(CTMDiagnostic(**d) for d in data)

    def save_diaginfo(self, filename):
        """
        TODO:
        """
        raise NotYetImplementedError()

    def save_tracerinfo(self, filename):
        """
        TODO:
        """
        raise NotYetImplementedError()

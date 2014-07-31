# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013-2014 Benoît Bovy
# see license.txt for more details
#

"""
Read / write 'diaginfo.dat' and 'tracerinfo.dat' files (GEOS-Chem diagnostics
metadata).

"""

from pygchem.utils.tff import read_fmt_file, write_fmt_file

DIAGINFO_FMT = (
    ('offset', 0, 8, int),
    ('name', 9, 40, str),
    ('description', 49, 9999, str)
)
TRACERINFO_FMT = (
    ('name', 0, 8, str),
    ('full_name', 9, 30, str),
    ('molecular_weight', 39, 10, float),
    ('carbon_weight', 49, 3, int),
    ('number', 52, 9, int),
    ('scale', 61, 10, float),
    ('unit', 72, 40, str)
)


def read_diaginfo(filename, **kwargs):
    """
    Read a 'diaginfo.dat' file and Returns a dictionary.
    """
    return read_fmt_file(filename, DIAGINFO_FMT, **kwargs)


def read_tracerinfo(filename, **kwargs):
    """
    Read a 'tracerinfo.dat' file and Returns a dictionary.
    """
    return read_fmt_file(filename, TRACERINFO_FMT, **kwargs)


def write_diaginfo(diag_metadata, filename):
    """
    Write diagnostics metadata to a 'diaginfo.dat' file.
    """
    return write_fmt_file(diag_metadata, filename, DIAGINFO_FMT)


def write_tracerinfo(tracer_metadata, filename):
    """
    Write diagnostics metadata to a 'diaginfo.dat' file.
    """
    return write_fmt_file(tracer_metadata, filename, TRACERINFO_FMT)

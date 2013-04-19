# -*- coding: utf-8 -*-

# module io.bpch2
# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2012-2013 Gerrit Kuhlmann, Benoît Bovy
# see license.txt for more details
# 
#
# Last modification: 04/2013

"""
Module for reading and writing binary punch files v2
"""

import os

import numpy as np

from pygchem.utils import uff, timetools


def read_bpch2(filename, mode='rb', endian='>', skip_values=True):
    """
    Read the binary punch file v2 given by `filename`, with `mode` and
    `endian`.
    
    Return ctm_file (the opened file instance), filetype, filetitle and
    the list of data blocks (both metadata and values as a dictionary).
    
    If `skip_values' is True, only data block headers will be read (values not
    loaded into memory). 
    """

    datablocks = []

    ctm_file = uff.FortranFile(filename, mode, endian)  # don't close file

    filetype = ctm_file.readline().strip()
    fsize = os.path.getsize(filename)
    filetitle = ctm_file.readline().strip()

    while ctm_file.tell() < fsize:
        # read first and second header line
        line = ctm_file.readline('20sffii')
        modelname, res0, res1, halfpolar, center180 = line
        line = ctm_file.readline('40si40sdd40s7i')
        category, index, unit, tau0, tau1, reserved = line[:6]
        dim0, dim1, dim2, dim3, dim4, dim5, skip = line[6:]

        # skip datablock (or read datablock if accessing value)
        position = ctm_file.tell()
        if skip_values:
            ctm_file.skipline()
            values = np.array([])
        else:
            values = np.array(ctm_file.readline('*f'))
            values = values.reshape((dim0, dim1, dim2), order='F')

        datablock = {'index': int(index),
                     'category': category.strip(),
                     'times': (timetools.tau2time(tau0),
                               timetools.tau2time(tau1)),
                     'modelname': modelname.strip(),
                     'center180': bool(center180),
                     'halfpolar': bool(halfpolar),
                     'origin': (dim3, dim4, dim5),
                     'resolution': (res0, res1),
                     'shape': (dim0, dim1, dim2),
                     'ctm_file': ctm_file,
                     'position': position,
                     'values': values,
                     'unit': unit.strip()
                    }

        datablocks.append(datablock)

    return ctm_file, filetype, filetitle, datablocks


def write_bpch2(ctm_file, filename, endian='>'):
    """
    Save a CTM file object to a binary punch file v2 given by `filename`. 
    """
    with uff.FortranFile(filename, 'wb', endian) as out_file:
        # write header
        out_file.writeline('40s', ctm_file.filetype.ljust(40))
        out_file.writeline('80s', ctm_file.title.ljust(80))

        # write data blocks
        for b in ctm_file.datablocks:
            out_file.writeline('20sffii',
                b.modelname.ljust(20), b.resolution[0],
                b.resolution[1], b.halfpolar, b.center180
                )
            out_file.writeline('40si40s2d40s7i',
                b.category.ljust(40), b.index, b.unit.ljust(40),
                timetools.time2tau(b.times[0]), timetools.time2tau(b.times[1]),
                ''.ljust(40), b.shape[0], b.shape[1], b.shape[2],
                b.origin[0], b.origin[1], b.origin[2], b.size * 4
                )
            vals = b.values.flatten('F')
            out_file.writeline('%df' % b.size, *vals)


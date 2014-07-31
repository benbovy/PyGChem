# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2012-2014 Gerrit Kuhlmann, Benoît Bovy
# see license.txt for more details
#

"""
Read / write binary punch (BPCH) files.

"""

import os

import numpy as np

from pygchem.utils import uff, timetools


FILETYPE02 = "CTM bin 02"
FILETYPE4D = "CTM bin 4D"
DEFAULT_TITLE = "GEOS-CHEM binary punch file v. 2.0"


def read_bpch(filename, mode='rb', skip_data=True, **kwargs):
    """
    Read the binary punch file v2 format.

    Parameters
    ----------
    filename : string
        name or path to the bpch file.
    mode : {'r', 'r+', rb', 'r+b', 'a'}
        file open mode (see :func:`open`). Writing only ('w' or 'wb') is not
        allowed.
    skip_data : bool
        if True, only data block metadata will be read (it will not load data
        into memory but data position and size information will be provided
        in addition to metadata).
    **kwargs
        extra parameters passed to :class:`pygchem.utils.uff.FortranFile`
        (e.g., `endian`).
    
    Returns
    -------
    bpch_file
        the open file instance (:class:`pygchem.utils.uff.FortranFile` object).
    filetype
        bpch file type identifier (given in the file's header)
    filetitle
        title (given in the file's header)
    datablocks
        the list of data blocks, i.e., dictionaries with data block metadata
        (and data).

    """
    if mode != 'a' and not mode.endswith('b'):
        mode += 'b'      # platform independent
    if 'w' in mode:
        raise ValueError("write-only mode is not allowed for reading the "
                         "bpch file")
    bpch_file = uff.FortranFile(filename, mode, **kwargs)  # don't close file

    datablocks = []
    filetype = bpch_file.readline().strip()
    fsize = os.path.getsize(filename)
    filetitle = bpch_file.readline().strip()

    while bpch_file.tell() < fsize:
        # read first and second header line
        line = bpch_file.readline('20sffii')
        modelname, res0, res1, halfpolar, center180 = line
        line = bpch_file.readline('40si40sdd40s7i')
        category, number, unit, tau0, tau1, reserved = line[:6]
        dim0, dim1, dim2, dim3, dim4, dim5, skip = line[6:]

        # skip datablock (or read datablock if accessing data)
        file_position = bpch_file.tell()
        if skip_data:
            bpch_file.skipline()
            data = np.array([])
        else:
            data = np.array(bpch_file.readline('*f'))
            data = data.reshape((dim0, dim1, dim2), order='F')

        datablock = {'number': int(number),
                     'category': category.strip(),
                     'times': (timetools.tau2time(tau0),
                               timetools.tau2time(tau1)),
                     'modelname': modelname.strip(),
                     'center180': bool(center180),
                     'halfpolar': bool(halfpolar),
                     'origin': (dim3, dim4, dim5),
                     'resolution': (res0, res1),
                     'shape': (dim0, dim1, dim2),
                     'from_file': os.path.abspath(filename),
                     'file_position': file_position,
                     'data': data,
                     'unit': unit.strip()}
        datablocks.append(datablock)

    return bpch_file, filetype, filetitle, datablocks


def create_bpch(filename, title=DEFAULT_TITLE, filetype=FILETYPE02, **kwargs):
    """
    Create a new empty bpch file.

    Parameters
    ----------
    filename : string
        name or path to the bpch file.
    title : string
        a title line to write in the file's header.
    filetype : string
        bpch file type identifier (either :attr:`bpch.FILETYPE02` or
        :attr:`bpch.FILETYPE4D`).
    **kwargs
        extra parameters passed to :class:`pygchem.utils.uff.FortranFile`
        (e.g., `endian`).

    Returns
    -------
    bpch_file
        the open file instance (:class:`pygchem.utils.uff.FortranFile` object).

    """
    bpch_file = uff.FortranFile(filename, 'wb', **kwargs)
    bpch_file.writeline('40s', filetype.ljust(40))
    bpch_file.writeline('80s', title.ljust(80))
    return bpch_file


def append_bpch(bpch_file, datablock):
    """
    Append a data block to an open bpch file.

    Parameters
    ----------
    bpch_file : file object
        bpch file (with writing permissions), as returned by
        :func:`read_bpch` or :func:`create_bpch`.
    datablock : dict
        data block metadata and data.

    """
    bpch_file.writeline(
        '20sffii',
        datablock['modelname'].ljust(20),
        datablock['resolution'][0], datablock['resolution'][1],
        datablock['halfpolar'], datablock['center180']
    )
    bpch_file.writeline(
        '40si40s2d40s7i',
        datablock['category'].ljust(40),
        datablock['number'], datablock['unit'].ljust(40),
        timetools.time2tau(datablock['times'][0]),
        timetools.time2tau(datablock['times'][1]),
        ''.ljust(40),
        datablock['shape'][0], datablock['shape'][1], datablock['shape'][2],
        datablock['origin'][0], datablock['origin'][1], datablock['origin'][2],
        datablock['data'].size * 4
    )
    data = datablock['data'].flatten('F')
    bpch_file.writeline('%df' % datablock['data'].size, *data)


def write_bpch(filename, datablocks, title=DEFAULT_TITLE,
               filetype=FILETYPE02, **kwargs):
    """
    Write data blocks to the binary punch file v2 format.

    Parameters
    ----------
    filename : string
        name or path to the bpch file.
    datablocks : sequence of dicts
        data blocks metadata and data.
    title : string
        a title line to write in the file's header.
    filetype : string
        bpch file type identifier (either :attr:`bpch.FILETYPE02` or
        :attr:`bpch.FILETYPE4D`).
    **kwargs
        extra parameters passed to :class:`pygchem.utils.uff.FortranFile`
        (e.g., `endian`).

    """
    with create_bpch(filename, title, filetype, **kwargs) as bpch_file:
        for db in datablocks:
            append_bpch(bpch_file, db)

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

from pygchem.datafield_backends.diagnostics import CTMDiagnosticInfo
from pygchem.utils import uff, timeutil, exceptions


FILETYPE02 = "CTM bin 02"
FILETYPE4D = "CTM bin 4D"
DEFAULT_TITLE = "GEOS-CHEM binary punch file v. 2.0"


class BPCHDataProxy(object):
    """A reference to the data payload of a single BPCH file datablock."""

    __slots__ = ('shape', 'dtype', 'path', 'endian', 'file_position',
                 'scale_factor', 'fill_value', '_data')

    def __init__(self, shape, dtype, path, endian, file_position,
                 scale_factor, fill_value):
        self.shape = shape
        self.dtype = dtype
        self.path = path
        self.fill_value = fill_value
        self.endian = endian
        self.file_position = file_position
        self.scale_factor = scale_factor
        self._data = None

    @property
    def ndim(self):
        return len(self.shape)

    @property
    def data(self):
        if self._data is None:
            self._data = self.load()
        return self._data

    def load(self):
        with uff.FortranFile(self.path, 'rb', self.endian) as bpch_file:
            bpch_file.seek(self.file_position)
            data = np.array(bpch_file.readline('*f'))
            data = data.reshape(self.shape, order='F')
        return data * self.scale_factor

    def __getitem__(self, keys):
        return self.data[keys]

    def __repr__(self):
        fmt = '<{self.__class__.__name__} shape={self.shape}' \
              ' dtype={self.dtype!r} path={self.path!r}' \
              ' file_position={self.file_position}' \
              ' scale_factor={self.scale_factor}>'
        return fmt.format(self=self)

    def __getstate__(self):
        return {attr: getattr(self, attr) for attr in self.__slots__}

    def __setstate__(self, state):
        for key, value in state.iteritems():
            setattr(self, key, value)


def read_bpch(filename, mode='rb', skip_data=True,
              diaginfo_file='', tracerinfo_file='', **kwargs):
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
        in addition to metadata, so that data can further be easily loaded).
    diaginfo_file : string
        path to the 'diaginfo.dat' file (optional). If empty, it will look
        for 'diaginfo.dat' in the same directory than `filename` or it will
        take a default one.
    tracerinfo_file : string
        path to the 'tracerinfo.dat' file (or empty string).
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
        and data (or a :class:`BPCHDataProxy` instance if `skip_data` is True).

    """
    if mode != 'a' and not mode.endswith('b'):
        mode += 'b'      # platform independent
    if 'w' in mode:
        raise ValueError("write-only mode is not allowed for reading the "
                         "bpch file")

    dir_path = os.path.dirname(filename)
    if not dir_path:
        dir_path = os.getcwd()
    if not tracerinfo_file:
        tracerinfo_file = os.path.join(dir_path, "tracerinfo.dat")
        if not os.path.exists(tracerinfo_file):
            tracerinfo_file = ''
    if not diaginfo_file:
        diaginfo_file = os.path.join(dir_path, "diaginfo.dat")
        if not os.path.exists(diaginfo_file):
            diaginfo_file = ''
    ctm_info = CTMDiagnosticInfo(diaginfo_file=diaginfo_file,
                                 tracerinfo_file=tracerinfo_file)

    with uff.FortranFile(filename, mode, **kwargs) as bpch_file:
        datablocks = []
        filetype = bpch_file.readline().strip()
        fsize = os.path.getsize(filename)
        filetitle = bpch_file.readline().strip()

        while bpch_file.tell() < fsize:
            # read first and second header line
            line = bpch_file.readline('20sffii')
            modelname, res0, res1, halfpolar, center180 = line
            line = bpch_file.readline('40si40sdd40s7i')
            category_name, number, unit, tau0, tau1, reserved = line[:6]
            dim0, dim1, dim2, dim3, dim4, dim5, skip = line[6:]

            # get additional metadata from tracerinfo / diaginfo
            try:
                cat = ctm_info.categories.select_item(category_name.strip())
                cat_attr = cat.to_dict()
                diag = ctm_info.diagnostics.select_item(
                    cat.offset + int(number)
                )
                diag_attr = diag.to_dict()
            except exceptions.SelectionMismatchError:
                diag = {'name': '', 'scale': 1}
                diag_attr = 'no additional metadata found for tracer/diagnostic'
                cat_attr = 'no metadata found for diagnostic category'

            # parse metadata, get data or set a data proxy
            if dim2 == 1:
                data_shape = (dim0, dim1)         # 2D field
            else:
                data_shape = (dim0, dim1, dim2)
            from_file = os.path.abspath(filename)
            file_position = bpch_file.tell()
            if skip_data:
                bpch_file.skipline()
                data = BPCHDataProxy(data_shape, 'f',
                                     from_file, bpch_file.endian,
                                     file_position, diag['scale'], np.nan)
            else:
                data = np.array(bpch_file.readline('*f'))
                data = data.reshape((dim0, dim1, dim2), order='F')

            datablock = {'number': int(number),
                         'name': diag['name'],
                         'category': category_name.strip(),
                         'times': (timeutil.tau2time(tau0),
                                   timeutil.tau2time(tau1)),
                         'modelname': modelname.strip(),
                         'center180': bool(center180),
                         'halfpolar': bool(halfpolar),
                         'origin': (dim3, dim4, dim5),
                         'resolution': (res0, res1),
                         'shape': data_shape,
                         'from_file': from_file,
                         'file_position': file_position,
                         'data': data,
                         'unit': unit.strip(),
                         'tracerinfo': diag_attr,
                         'diaginfo': cat_attr}
            datablocks.append(datablock)

    return filetype, filetitle, datablocks


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
    if isinstance(datablock['data'], BPCHDataProxy):
        data = datablock['data'].data / datablock['data'].scale_factor
    else:
        data = datablock['data']

    if len(datablock['shape']) == 2:
        data_shape = datablock['shape'] + (1,)    # 2D field
    else:
        data_shape = datablock['shape']

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
        timeutil.time2tau(datablock['times'][0]),
        timeutil.time2tau(datablock['times'][1]),
        ''.ljust(40),
        data_shape[0], data_shape[1], data_shape[2],
        datablock['origin'][0], datablock['origin'][1], datablock['origin'][2],
        data.size * 4
    )
    data_array = data.flatten('F')
    bpch_file.writeline('%df' % data.size, *data_array)


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

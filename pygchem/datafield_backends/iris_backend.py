# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013-2014 Benoit Bovy
# see license.txt for more details
#

"""
Iris backend (http://scitools.org.uk/iris).

Also adds support for BPCH format to Iris.

"""

from __future__ import absolute_import

import glob
from types import StringTypes

import iris
import iris.unit
import iris.coords
import iris.fileformats
import iris.io.format_picker as fp

from pygchem import grid
from pygchem.utils import uff
from pygchem.io import bpch
from pygchem.tools import irisutil


# -----------------------------------------------------------------------------
# BPCH Support for Iris
#------------------------------------------------------------------------------

def cubes_from_bpch(filenames, callback=None, **kwargs):
    """
    Return a generator of Iris cubes from BPCH filenames.

    Parameters
    ----------
    filenames : string or sequence of strings
        (list of) BPCH filename(s) to load (can be UNIX expressions,
        e.g., '*.bpch').
    callback : func
        a function which can be passed on to :func:`iris.io.run_callback`
    **kwargs
        any extra keyword argument passed to
        :class:`pygchem.utils.uff.FortranFile` (e.g., `endian`).

    Notes
    -----
    The resultant cubes may not be in the same order as in the files.

    """
    if isinstance(filenames, StringTypes):
        filenames = [filenames]

    for filename in filenames:
        for path in glob.glob(filename):
            
            filetype, filetitle, datablocks = bpch.read_bpch(path, **kwargs)
                
            # assume that CTM Grid is the same for all datablocks
            # (compute coordinates once per file to save CPU/memory).
            leading_datablock = datablocks[0]
            ctm_grid = grid.CTMGrid.from_model(
                leading_datablock['modelname'],
                resolution=leading_datablock['resolution']
            )
            lon, lat, lev = irisutil.coord_from_grid(
                ctm_grid, coord=('longitude', 'latitude', 'levels')
            )

            for datablock in datablocks:

                if len(datablock['shape']) == 2:
                    dim_coords = [(lon, 0), (lat, 1)]
                else:
                    dim_coords = [(lon, 0), (lat, 1), (lev, 2)]  # 3D default

                cube = irisutil._datablock_to_cube(
                    datablock,
                    dim_coords_and_dims=dim_coords,
                    coords_from_model=False,
                    errcoord='pass'
                )

                if callback is not None:
                    cube = iris.io.run_callback(callback,
                                                cube,
                                                datablock,
                                                filename)
                if cube is None:
                    continue
                yield cube

# TODO: custom loaders for ND49, ND50


class LeadingLineUff(fp.FileElement):
    """
    A :class:`FileElement` that returns the (unpacked) first line of an
    un-formatted binary Fortran file.

    """
    def __init__(self, endian):
        self._endian = endian
        super(LeadingLineUff, self).__init__(requires_fh=False)

    def get_element(self, basename, file_handle):
        with uff.FortranFile(basename, endian=self._endian) as uffile:
            leading_line = uffile.readline()
        return leading_line


bpch_spec = fp.FormatSpecification(
    'Binary Punch File (BPCH) v2 (big-endian)',
    LeadingLineUff('>'),
    lambda line: line.lstrip().startswith(bpch.FILETYPE02),
    cubes_from_bpch,
    priority=7
)
iris.fileformats.FORMAT_AGENT.add_spec(bpch_spec)
bpch_spec_le = fp.FormatSpecification(
    'Binary Punch File (BPCH) v2 (little-endian)',
    LeadingLineUff('<'),
    lambda line: line.lstrip().startswith(bpch.FILETYPE02),
    lambda filenames, callback=None: cubes_from_bpch(filenames, callback,
                                                     endian='<'),
    priority=6
)
iris.fileformats.FORMAT_AGENT.add_spec(bpch_spec_le)
# TODO: CTM bin 4D  format


# -----------------------------------------------------------------------------
# Implementation of CTM datafield classes and functions
#------------------------------------------------------------------------------

load = iris.load
load_cube = iris.load_cube
load_cubes = iris.load_cubes
load_raw = iris.load_raw
#load_strict = iris.load_strict    # depreciated

save = iris.save

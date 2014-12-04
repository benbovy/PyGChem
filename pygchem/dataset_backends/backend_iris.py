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
import collections
import warnings
import datetime
from types import StringTypes

import numpy as np
import netCDF4 as nc
import iris
import iris.cube
import iris.unit
import iris.coords
import iris.analysis
import iris.exceptions
import iris.fileformats
import iris.std_names
import iris.io.format_picker as fp

from pygchem import grid, diagnostics
from pygchem.utils import uff
from pygchem.io import bpch
from pygchem.tools import irisutil, ctm2cf, timeutil


# -----------------------------------------------------------------------------
# BPCH Support for Iris
# -----------------------------------------------------------------------------

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
            if leading_datablock['origin'] != (1, 1, 1):
                imin = np.array(leading_datablock['origin']) - 1
                imax = imin + np.array(leading_datablock['shape'])
                region_box = zip(imin, imax)
            else:
                region_box = None
            lon, lat, lev = irisutil.coord_from_grid(ctm_grid,
                                                     region_box=region_box)

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

                # add `filetitle` in the cube's attributes
                # (commented: may not properly concatenate cubes)
                #cube.attributes['bpch_title'] = filetitle

                if callback is not None:
                    cube = iris.io.run_callback(callback,
                                                cube,
                                                datablock,
                                                filename)
                if cube is None:
                    continue
                yield cube


class LeadingLineUff(fp.FileElement):
    """
    A :class:`FileElement` child that returns the (unpacked) first line of an
    un-formatted binary Fortran file.

    """
    def __init__(self, endian):
        self._endian = endian
        super(LeadingLineUff, self).__init__(requires_fh=True)

    def get_element(self, basename, file_handle):
        with uff.FortranFile(file_handle.name, endian=self._endian) as uffile:
            try:
                leading_line = uffile.readline()
            except IOError:
                leading_line = ''
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
# TODO: CTM bin 4D  format (?)

# TODO: write a bpch Saver class (not a priority)


# -----------------------------------------------------------------------------
# Callbacks for loading cubes (various CF-related fixes)
# -----------------------------------------------------------------------------

# Init caches for Iris coordinate objects (save memory and CPU)
_coordcache = dict()
_coordcache2 = dict()


def fix_bpch2nc(cube, field, filename):
    """
    An Iris load callback for properly loading the NetCDF files
    created by BPCH2NC (GAMAP v2-12+).

    """
    global _coordcache

    # units
    units = field.unit.strip()
    try:
        cube.units = units
    except ValueError:
        # Try to get equivalent units compatible with udunits.
        # Store original unit as cube attribute
        conform_units = ctm2cf.get_cfcompliant_units(units)
        try:
            cube.units = conform_units
        except ValueError:
            warnings.warn("Invalid udunits2 '{0}'".format(units))
    cube.attributes["ctm_units"] = units

    # a hack for keeping cube's long_name but show var_name in cube summary
    iris.std_names.STD_NAMES[cube.var_name] = {'canonical_units': cube.units}
    cube.standard_name = cube.var_name

    # add spatial coordinates
    modelname = cube.attributes.get('Model')
    res = cube.attributes.get('Delta_Lon'), cube.attributes.get('Delta_Lat')
    nlayers = cube.attributes.get('NLayers')
    cache_key = modelname, res, nlayers

    if _coordcache.get(cache_key) is None:
        ctm_grid = grid.CTMGrid.from_model(modelname, resolution=res)
        coord_names = ('longitude', 'latitude', 'levels')
        lon, lat, lev = irisutil.coord_from_grid(ctm_grid, coord=coord_names)
        _coordcache[cache_key] = {'lon': lon, 'lat': lat, 'lev': lev}
    if cube.ndim == 3:
        cube.add_dim_coord(_coordcache[cache_key]['lon'], 2)
        cube.add_dim_coord(_coordcache[cache_key]['lat'], 1)
        if cube.shape[0] == nlayers:
            cube.add_dim_coord(_coordcache[cache_key]['lev'], 0)

    # add time scalar coordinate (get info from global attributes)
    tstart = (str(cube.attributes['Start_Date']) +
              str(cube.attributes['Start_Time']))
    tstart = datetime.datetime.strptime(tstart, "%Y%m%d%H")
    tstart = timeutil.time2tau(tstart)
    tend = (str(cube.attributes['End_Date']) +
            str(cube.attributes['End_Time']))
    tend = datetime.datetime.strptime(tend, "%Y%m%d%H")
    tend = timeutil.time2tau(tend)
    time_coord = iris.coords.DimCoord(points=tstart,
                                      bounds=(tstart, tend),
                                      standard_name='time',
                                      units=irisutil.CTM_TIME_UNIT_IRIS)
    cube.add_aux_coord(time_coord)

    # attributes
    # TODO: don't remove all attributes
    cube.attributes.clear()

    # handle 2D fields (remove 1-sized 1st dimension)
    if cube.ndim == 3 and cube.shape[0] == 1:
        #dummy_coord = iris.coords.DimCoord([0], long_name='dummy',
        #                                   bounds=[0, 1])
        #cube.add_dim_coord(dummy_coord, 0)
        return cube.slices(range(1, cube.ndim)).next()


def fix_bpch2coards(cube, field, filename):
    """
    An Iris load callback for properly loading the NetCDF files
    created by BPCH2COARDS (GAMAP v2-17+).

    """
    global _coordcache2

    # units
    units = field.units
    try:
        cube.units = units
    except ValueError:
        # Try to get equivalent units compatible with udunits.
        # Store original unit as cube attribute
        conform_units = ctm2cf.get_cfcompliant_units(units)
        try:
            cube.units = conform_units
        except ValueError:
            warnings.warn("Invalid udunits2 '{0}'".format(units))
    cube.attributes["ctm_units"] = units

    # a hack for keeping cube's long_name but show var_name in cube summary
    iris.std_names.STD_NAMES[cube.var_name] = {'canonical_units': cube.units}
    cube.standard_name = cube.var_name

    # attributes
    # TODO: don't remove all attributes
    cube.attributes.clear()

    # longitude coordinate (non strictly monotonic) degrees -> degrees_east
    try:
        lon = cube.coord('longitude')
        lon_dim = cube.coord_dims(lon)[0]
        cache_key = 'longitude', filename

        if _coordcache2.get(cache_key) is None:
            west_ind = np.nonzero(lon.points >= 180.)
            lon.points[west_ind] = -1. * (360. - lon.points)
            lon.units = 'degrees_east'
            _coordcache2[cache_key] = iris.coords.DimCoord.from_coord(lon)

        cube.remove_coord(lon)
        cube.add_dim_coord(_coordcache2[cache_key], lon_dim)
    except iris.exceptions.CoordinateNotFoundError:
        pass

    # levels coordinate
    # 'sigma_level' depreciated in the CF standard (not supported by UDUNITS)
    try:
        lev = cube.coord('Eta Centers')
        lev_dim = cube.coord_dims(lev)[0]
        lev_std_name = 'atmosphere_hybrid_sigma_pressure_coordinate'
        cache_key = lev_std_name, filename

        if _coordcache2.get(cache_key) is None:
            lev.standard_name = lev_std_name
            lev.units = iris.unit.Unit('1')
            d = nc.Dataset(filename)
            elev = d.variables['edge'][:]
            lev.bounds = np.column_stack((elev[:-1], elev[1:]))
            _coordcache2[cache_key] = iris.coords.DimCoord.from_coord(lev)

        cube.remove_coord(lev)
        cube.add_dim_coord(_coordcache2[cache_key], lev_dim)
    except iris.exceptions.CoordinateNotFoundError:
        pass

    # time: dimension -> scalar coordinate (+ add bounds)
    try:
        time_coord = cube.coord('time')
        time_dim = cube.coord_dims(time_coord)[0]

        with iris.FUTURE.context(cell_datetime_objects=True):
            tstart = time_coord.cell(0).point
        delta_t = time_coord.attributes.pop('delta_t')
        tend = tstart + timeutil.strp_relativedelta(delta_t)
        time_coord.bounds = [timeutil.time2tau(tstart),
                             timeutil.time2tau(tend)]
        if cube.shape[time_dim] == 1:
            slices_dims = [d for d in range(cube.ndim) if d != time_dim]
            return cube.slices(slices_dims).next()
    except iris.exceptions.CoordinateNotFoundError:
        pass


# -----------------------------------------------------------------------------
# Wrappers for Iris load and save functions
# -----------------------------------------------------------------------------

load = iris.load
load_dataset = iris.load_cube
load_datasets = iris.load_cubes
load_raw = iris.load_raw

load_callbacks = {
    'gamap_bpch2nc': fix_bpch2nc,
    'gamap_bpch2coards': fix_bpch2coards
}

save = iris.save

# TODO: an alternative way to link CTM metadata, other than cube attributes ?
#   - cube attributes may +- slow the cube creation and concatenation
#   - CTM metadata related to the geometry, useful because of the metadata
#     format but duplicates the information in the cube
#
#   - solution? add attributes 'diagnostic', 'bpch_metadata'
#               and overriding the copy method of the cube?

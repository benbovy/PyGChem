# -*- coding: utf-8 -*-

# module iris
# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013 Benoit Bovy
#
#
# see license.txt for more details
#
#
# Last modification: 08/2013

"""
Iris backend (http://scitools.org.uk/iris).

Also adds support for BPCH format to Iris.

"""

#TODO: make custom 'load', 'load_cubes' and 'load_cube' functions on top of iris
#equivalent funcs, which allow user to specify additional information. Useful
#for example in the case of ND49 and ND50 diagnostic files which can consist of
#CTM Grid subsets with unknown bounds.
#i.e. : load('file.bpch', **kwargs)
#---> to complicated: let user load ND49 and ND50 as datablocks (using the
# simple, low-level API) and then convert to cubes using datablocks_to_cubes
# with appropriate dim coordinates. (show this in examples).  

from __future__ import absolute_import

import os
import glob
import warnings
from copy import deepcopy

import numpy as np
import iris
import iris.exceptions as iexceptions
import iris.coords as icoords
import iris.fileformats as ifileformats
import iris.io.format_picker as format_picker

from pygchem.grid import CTMGrid
from temp.diagnostics import CTMFile, DataBlock
from pygchem.datafield_backends import gc2cf
from pygchem.utils import uff, timetools



# time encoding for bpch files
TIME_UNIT = iris.unit.Unit('hours since 1985-01-01 00:00:00',
                           calendar=iris.unit.CALENDAR_STANDARD)


class _BPCHDataProxy(object):
    """
    A reference to the data payload of a single field in a BPCH file.
    """

    __slots__ = ('path', 'variable_name', 'position', 'endian', 'scale_factor')

    def __init__(self, path, variable_name, position, endian, scale_factor):
        self.path = path
        self.variable_name = variable_name
        self.position = position
        self.endian = endian
        self.scale_factor = scale_factor

    def __repr__(self):
        return '%s(%r, %r at %d)' % (self.__class__.__name__, self.path,
                               self.variable_name, self.position)

    def __getstate__(self):
        return {attr: getattr(self, attr) for attr in self.__slots__}

    def __setstate__(self, state):
        for key, value in state.iteritems():
            setattr(self, key, value)

    def load(self, data_shape, data_type, mdi, deferred_slice):
        """
        Load the corresponding proxy data item and perform any deferred
        slicing.

        Parameters
        ----------
        data_shape : tuple of int
            The data shape of the proxy data item.
        data_type : :class:`numpy.dtype`
            The data type of the proxy data item.
        mdi : float
            The missing data indicator value.
        deferred_slice : tuple
            The deferred slice to be applied to the proxy data item.

        Returns
        -------
            :class:`numpy.ndarray`

        """
        with uff.FortranFile(self.path, 'rb', self.endian) as bpch_file:
            bpch_file.seek(self.position)
            variable = np.array(bpch_file.readline('*f'))
            variable = variable.reshape(data_shape, order='F')
            payload = variable[deferred_slice] * self.scale_factor

        return payload


def _cube_data_reshape(cube, new_shape):
    """
    Reshape the data of `cube`, given a `new_shape`.
    
    This function is only applicable for adding or removing one-size dimensions
    to the data.
    
    See :func:`remove_one_size_dims` and :func:`add_one_size_dims`.
    """
    
    if cube._data_manager is not None:
        
        # create a new iris.fileformats.manager.DataManager instance.
        # and assign it to the cube
        # `new_shape` is the cube's shape = proxy_array shape + data
        # shape --> retreive the 'data member' of the new shape.
        
        data_manager = cube._data_manager
        
        old_shape = data_manager._orig_data_shape
        new_shape = new_shape[len(cube._data.shape):]
        
        if np.prod(old_shape) != np.prod(new_shape):
            raise ValueError("Given new shape is not compatible with "
                             "array shape in file")
        
        new_data_manager = ifileformats.manager.DataManager(
                              data_shape=new_shape,
                              data_type=deepcopy(data_manager.data_type),
                              mdi=deepcopy(data_manager.mdi))
        cube._data_manager = new_data_manager
    
    else:
        cube.data = cube.data.reshape(new_shape)


def remove_one_size_dims(cube):
    """
    Remove all one-size dimensions of the cube (in-place).

    Parameters
    ----------
    cube : :class:`iris.cube.Cube` instance
        the Iris's that will be modified in place.

    Notes
    -----
    Any dimension coordinate mapped to a 1-size dimension
    will be re-defined as an auxiliary coordinate.

    No copy nor access to data is made. 
    """
    
    new_shape = tuple((n for n in cube.shape if n > 1))
    
    remove_dims = [True if n == 1 else False for n in cube.shape]
    coord_dims = [None] * cube.ndim
    
    for coord, dim in cube._dim_coords_and_dims:
        coord_dims[dim] = coord
    
    for remove_dim, coord_dim in zip(remove_dims, coord_dims):
        if remove_dim:
            if coord_dim is not None:
                cube.remove_coord(coord_dim)
                cube.add_aux_coord(coord_dim)
    
    _cube_data_reshape(cube, new_shape)


def add_one_size_dims(cube, dim_coords_and_dims):
    """
    Not implemented yet !!
    
    Add one-size dimension(s) to the cube (in-place).
    
    """
    # TODO: 
    raise iexceptions.NotYetImplementedError("this function is not "
                                             "implemented yet")


def levels_as_aux_factories(grid_box_height, bottom_lev_pressure,
                            orography=None):
    """
    Return auxiliary coordinate factories for vertical levels from pressure
    and grid box height variables given as Iris cubes.
    
    Parameters
    ----------
    grid_box_height : :class:`iris.cube.Cube` instance
        Height of the grid boxes. Does correspond to the
        **BXHEIGHT** diagnostic of the **BXHGHT-$** category.
    bottom_lev_pressure : :class:`iris.cube.Cube` instance
        Pressure at the bottom edge of the grid boxes. Does correspond to the
        **PSURF** diagnostic of the **PEDGE-$** category.
    orography : :class:`iris.cube.Cube` instance or None
        Orography (height of the Geoid) (optional).
    
    Returns
    -------
    A list of :class:`iris.aux_factory.AuxCoordFactory` subclasses which
    can be used to build auxiliary vertical coordinates on demand (e.g.,
    atmosphere_sigma_coordinate, atmosphere_hybrid_height_coordinate,
    atmosphere_hybrid_sigma_pressure_coordinate...).
    
    Only 'atmosphere_hybrid_sigma_pressure_coordinate' is currently
    implemented !
        
    Notes
    -----
    Obviously, all input cubes must be compatible.
    
    `bottom_lev_pressure` is used to get only surface pressures.
    
    """
    
    height_model = grid_box_height.attributes.get('modelname')
    press_model = bottom_lev_pressure.attributes.get('modelname')
    
    if height_model != press_model:
        raise iexceptions.InvalidCubeError("grid_box_height and "
                                           "bottom_lev_pressure must "
                                           "correspond to the same model")
    
    else:
        ctm_grid = CTMGrid.from_model(height_model)
    
    # Hybrid Pressure Factory
    # TODO:


def coords_from_ctm_grid(ctm_grid, **kwargs):
    """
    Get Iris 1-dimensional coordinates objects from a CTM grid.
    
    Parameters
    ----------
    ctm_grid : string or :class:`pygchem.grid.CTMGrid` instance
        Name of the CTM grid (model) or CTMGrid object.
    
    Returns
    -------
    lon_coord : :class:`iris.coords.DimCoord` instance
        Longitude coordinates
    lat_coord : :class:`iris.coords.DimCoord` instance
        Latitude coordinates
    lev_coord : :class:`iris.coords.DimCoord` instance or None 
        Vertical level coordinates (layer numbers), or None if `ctm_grid` has
        no vertical layers.
    sigma_coord : :class:`iris.coords.DimCoord` instance or None 
        Vertical level coordinates (sigma), or None if `ctm_grid` has
        no vertical layers or is hybrid.
    eta_coord : :class:`iris.coords.DimCoord` instance or None 
        Vertical level coordinates (ETA), or None if `ctm_grid` has
        no vertical layers or is not hybrid.
    press_coord : :class:`iris.coords.DimCoord` instance or None
        Vertical pressure levels, or None if `ctm_grid` has
        no vertical layers.
    height_coord : :class:`iris.coords.DimCoord` instance or None
        Vertical altitude levels, or None if `ctm_grid` has
        no vertical layers.
    
    Other Parameters
    ----------------
    **kwargs
        May be optional parameter(s) for :func:`grid.CTMGrid.from_model`.
    
    Notes
    -----
    The returned `sigma_coord`, `eta_coord`, `press_coord` and `height_coord`
    are only approximations of the real coordinates values, given a reference
    atmosphere surface pressure and a default model of pressure vs. altitude
    for `height_coord` (see :meth:`grid.CTMGrid.get_layers`).
    """
    
    if not isinstance(ctm_grid, CTMGrid):
        ctm_grid = CTMGrid.from_model(ctm_grid, **kwargs)
    
    def reshape_bounds(bounds_array):
        """Reshape bounds from shape (Nlayers + 1) to shape (Nlayers, 2)"""
        return np.column_stack((bounds_array[:-1], bounds_array[1:]))
    
    clon, clat = ctm_grid.lonlat_centers
    elon, elat = ctm_grid.lonlat_edges
    elon = reshape_bounds(elon)
    elat = reshape_bounds(elat)
    
    lon_coord = icoords.DimCoord(clon, standard_name="longitude",
                                 units="degrees_east", bounds=elon)
    lat_coord = icoords.DimCoord(clat, standard_name="latitude",
                                 units="degrees_north", bounds=elat)
    
    if ctm_grid.Nlayers is not None:
        levels = np.arange(1, ctm_grid.Nlayers + 1)
        lev_coord = icoords.DimCoord(levels,
                                     standard_name="model_level_number",
                                     units="1")
        
        if ctm_grid.hybrid:
            ceta = ctm_grid.eta_centers
            eeta = ctm_grid.eta_edges
            eeta = reshape_bounds(eeta)
            eta_coord = icoords.DimCoord(
                            ceta,
                            var_name="atmosphere_eta_coordinate_approx",
                            units="1",
                            bounds=eeta)
            sigma_coord = None
        else:
            csigma = ctm_grid.sigma_centers
            esigma = ctm_grid.sigma_edges
            esigma = reshape_bounds(esigma)
            sigma_coord = icoords.DimCoord(
                            csigma,
                            var_name="atmosphere_sigma_coordinate_approx",
                            units="1",
                            bounds=esigma)
            eta_coord = None
        
        cpress = ctm_grid.pressure_centers
        epress = ctm_grid.pressure_edges
        epress = reshape_bounds(epress)
        press_coord = icoords.DimCoord(
                cpress,
                var_name="atmosphere_hybrid_sigma_pressure_coordinate_approx",
                units="hPa",
                bounds=epress)
        
        cheight = ctm_grid.altitude_centers
        eheight = ctm_grid.altitude_edges
        eheight = reshape_bounds(eheight)
        height_coord = icoords.DimCoord(
                cheight,
                var_name="atmosphere_hybrid_height_coordinate_approx",
                units="km",
                bounds=eheight)
        
    else:
        lev_coord = None
        sigma_coord = None
        eta_coord = None
        press_coord = None
        height_coord = None
    
    return (lon_coord, lat_coord, lev_coord,
            sigma_coord, eta_coord,
            press_coord, height_coord)


def _datablock_to_cube(datablock, dim_coords_and_dims=None,
                       aux_coords_and_dims=None, aux_factories=None,
                       coords_from_model=True, error_coord=True,
                       reduce_dims=True, **kwargs):
    """
    Create a :class:`iris.cubes.Cube` instance given a
    :class:`pygchem.diagnostics.DataBlock` instance.
    
    See docstrings of :func:`datablocks_to_cubes` for explanations about the
    arguments.
    """
    
    # name to assign to the cube
    cf_name = "_".join([datablock.name, datablock.category])
    
    # Create cube with data proxy if datablock was loaded from a file...
    bpch_file = datablock._ctm_file
    
    if bpch_file:     
        filepath = os.path.abspath(bpch_file.name)
         
        data_proxy = np.array(_BPCHDataProxy(filepath,
                                             cf_name,
                                             datablock._position,
                                             bpch_file.endian,
                                             datablock.scale))
          
        dummy_data = np.zeros(1, dtype='f')            # Figure out what the
        if hasattr(datablock, 'scale'):                # data type will be
            dummy_data = datablock.scale * dummy_data  # after scale transform
         
        data_manager = ifileformats.manager.DataManager(datablock.shape,
                                                        dummy_data.dtype,
                                                        None)
         
        cube = iris.cube.Cube(data_proxy, data_manager=data_manager)
    
    # ... or create cube with data otherwise.
    else:
        cube = iris.cube.Cube(datablock.values)
    
    # set the cube's name
    cube.rename(cf_name)
    
    # set coordinates
    def add_coord(coord, dim, add_coord_func, error_coord=True):
        """add dim or aux coord to the cube (conditional)"""
        try:
            if coord is not None:
                add_coord_func(coord, dim)
        except ValueError:
            if error_coord:
                raise
            else:
                pass    # TODO: display warn or not ??
     
    # set coordinates from datablock metadata, using class grid.CTMGird
    if coords_from_model:
        ctm_grid = CTMGrid.from_model(datablock.modelname,
                                      resolution=datablock.resolution)
        lon, lat, lev, sigma, eta, press, alt = coords_from_ctm_grid(ctm_grid)
        cube.add_dim_coord(lon, 0)
        cube.add_dim_coord(lat, 1)
        if cube.ndim > 2:
            add_coord(lev, 2, cube.add_dim_coord, error_coord=False)
            for coord in [sigma, eta, press, alt]:
                add_coord(coord, 2, cube.add_aux_coord,
                          error_coord=error_coord)
     
    # add given dimension coordinates if any
    if dim_coords_and_dims:
        for coord, dim in dim_coords_and_dims:
            add_coord(coord, dim, cube.add_dim_coord, error_coord=error_coord)
     
    # add given auxiliary coordinates if any
    if aux_coords_and_dims:
        for coord, dim in aux_coords_and_dims:
            add_coord(coord, dim, cube.add_aux_coord, error_coord=error_coord)
     
    # add given aux coordinates factories if any
    if aux_factories:
        for aux_factory in aux_factories:
            cube.add_aux_factory(aux_factory)
     
    # time coordinates
    point = timetools.time2tau(datablock.times[0])
    bounds = [timetools.time2tau(time) for time in datablock.times]
 
    time_coord = icoords.DimCoord(points=point,
                                  bounds=bounds,
                                  standard_name='time',
                                  units=TIME_UNIT)
 
    cube.add_aux_coord(time_coord)
     
    # units
    units = datablock.unit.strip()
    try:
        cube.units = units
    except ValueError:
        # Try to get equivalent units compatible with udunits.
        # Store original unit as cube attribute
        conform_units = gc2cf.get_conforming_units(units)
        try:
            cube.units = conform_units
        except ValueError:
            warnings.warn("Unhandled units '{0}' recorded in cube attributes.".
                          format(units))
        cube.attributes["geoschem_units"] = units
     
    # attributes: Add datablock's attributes to the cube
    attrs = ['name', 'full_name', 'index', 'number', 'category', 'modelname',
             'chemical', 'molecular_weight', 'carbon_weight', 'hydrocarbon']
    for attr in attrs:
        if hasattr(datablock, attr):
            value = getattr(datablock, attr)
            cube.attributes[attr] = value
 
    # Copy datablock attributes that correspond to CF-reserved attributes
    # (e.g., datablock.scale -> CF 'scale_factor')
    # Iris doesn't (allow to) store these attributes in cube's metadata (a few
    # exceptions: history, units, standard_name...), but it can be useful in
    # some cases (e.g., I've read GEOS-Chem adjoint still works with
    # compressed/unscaled datablocks stored in input files ?).
    # workaround: store these attributes as a dictionary that is accessed by
    # defining 'cf_attrs' in cube's attributes.
    cf_attributes = dict()
    cf_attributes['scale_factor'] = datablock.scale 
    cube.attributes['cf_attrs'] = cf_attributes
    
    # add kwargs as attributes
    cube.attributes.update(kwargs)
    
    # remove one-size dimensions if desired
    if reduce_dims:
        remove_one_size_dims(cube)
    
    return cube


def datablocks_to_cubes(datablocks, dim_coords_and_dims=None,
                        aux_coords_and_dims=None, aux_factories=None,
                        coords_from_model=True, error_coord=True,
                        reduce_dims=True, merge=True, **kwargs):
    """
    Create Iris's cubes from BPCH datablocks.
    
    Parameters
    ----------
    datablocks : (sequence of) :class:`pygchem.diagnostics.DataBlock` instances
        Datablock(s) to convert to cubes
    dim_coords_and_dims : list or None
        A list of dimension coordinates with scalar dimension mappings, e.g
        ``[(lat_coord, 0), (lon_coord, 1)]``, relative to all datablocks.
        Coordinates must be :class:`iris.coords.DimCoord` instances.
    aux_coords_and_dims : list or None
        A list of auxiliary coordinates with scalar dimension mappings, e.g
        ``[(lat_coord2, 0), (lon_coord2, 1)]``, relative to all datablocks.
        Coordinates must be :class:`iris.coords.DimCoord` and/or
        :class:`iris.coords.AuxCoord` instances.
    aux_factories : list or None
        A list of auxiliary coordinate factories (see :mod:`iris.aux_factory`)
        relative to all datablocks.
    coords_from_model : bool
        If True, attributes of each datablock (e.g., `modelname`, `resolution`)
        will be used to set the cube coordinates
        (see :class:`pygchem.grid.CTMGrid` and :func:`coords_from_ctmgrid`).
    error_coord : bool
        If True, raise a ValueError when a given coordinate could not be mapped
        to the given dimension(s) of the cube (e.g., invalid or unavailable
        dimension).
        Otherwise, invalid (coordinate, dimension mapping) pairs are ignored.
    reduce_dims : bool
        If True, remove all unnecessary one-size dimensions
        (see :func:`remove_one_size_dims`).  
    merge : bool
        If True, contiguous datablocks in space or time are merged into a
        single cube when possible.
    
    Returns
    -------
    A :class:`iris.cubes.Cube` instance or a :class:`iris.cubes.CubeList`
    instance, following the given input.
    
    Other parameters
    ----------------
    **kwargs
        Any 'attribute_name'='attribute_value' that will be added as attributes
        for each returned cube.
    
    Notes
    -----
    The resultant cubes may not be in the same order as in the datablocks
    sequence.
    
    The names of the resultant cubes are given by "`name`_`category`"
    attributes of the corresponding datablocks.
    
    Examples
    --------
    TODO:
    """
    
    if isinstance(datablocks, DataBlock):
        datablocks = [datablocks]
    
    cube_gen = (_datablock_to_cube(datablock,
                                   dim_coords_and_dims,
                                   aux_coords_and_dims,
                                   aux_factories,
                                   coords_from_model,
                                   error_coord,
                                   reduce_dims,
                                   **kwargs)
                for datablock in datablocks)
    cube_list = iris.cube.CubeList(cube_gen)
    
    if merge:
        cube_list = cube_list.merge(unique=False)
    
    if len(cube_list) == 1:
        return cube_list[0]
    else:
        return cube_list
    

def _bpch_to_cubes(filenames, callback=None, endian='>'):
    """
    Return a generator of Iris cubes from BPCH filenames.

    Parameters
    ----------
    filenames : string or sequence of strings
        (list of) BPCH filename(s) to load (can be UNIX expressions,
        e.g., '*.bpch').
    callback : func
        a function which can be passed on to :func:`iris.io.run_callback`
    endian : {'@', '>', '<'}
        byte order of the binary files

    Notes
    -----
    The resultant cubes may not be in the same order as in the files.
    """
    
    if isinstance(filenames, basestring):
        filenames = [filenames]

    for filename in filenames:
        for path in glob.glob(filename):
            
            with CTMFile.fromfile(filename, mode='rb', endian=endian) as \
                       bpch_file:
                
                # to save memory and CPU, assume that CTM Grid is the
                # same for all datablocks in a given file. Compute coordinates
                # once per file.
                ref_datablock = bpch_file.datablocks[0]
                ctm_grid = CTMGrid.from_model(
                                        ref_datablock.modelname,
                                        resolution=ref_datablock.resolution)
                coords = coords_from_ctm_grid(ctm_grid)
                lon, lat, lev, sigma, eta, press, alt = coords
                dim_coords = [(lon, 0), (lat, 1), (lev, 2)]
                aux_coords = [(sigma, 2), (eta, 2), (press, 2), (alt, 2)] 
                
                for datablock in bpch_file.datablocks:
                
                    cube = _datablock_to_cube(datablock,
                                              dim_coords_and_dims=dim_coords,
                                              aux_coords_and_dims=aux_coords,
                                              coord_from_model=False,
                                              error_coord=False)

                    if callback is not None:
                        cube = iris.io.run_callback(callback,
                                                    cube,
                                                    datablock,
                                                    filename)
                    if cube is None:
                        continue
                    yield cube


def _bpch_to_cubes_le(filenames, callback=None):
    """
    Return a generator of Iris cubes from BPCH filenames
    (little endian only).
    
    See :func:`_bpch_as_cubes`.
    """ 
    return _bpch_to_cubes(filenames, callback, endian='<')


# Add BPCH file format (big endian, little endian) to auto-detection of file
# formats implemented in Iris, giving it a prior detection before other formats
def _read_bin_header(filename, endian):
    """
    Return the header (i.e., 1st line) of an unformatted binary Fortran file.
    """
    with uff.FortranFile(filename, endian=endian) as ufile:
        header = ufile.readline()
    return header

_LEADING_LINE_BIN = format_picker.FileElement(
                        'Leading line Binary',
                        lambda filename, fh: _read_bin_header(filename, '>'))
_LEADING_LINE_BIN_LE = format_picker.FileElement(
                        'Leading line Binary (little endian)',
                        lambda filename, fh: _read_bin_header(filename, '<'))

_BPCH_spec = format_picker.FormatSpecification(
                        'Binary Punch File (BPCH) v2',
                        _LEADING_LINE_BIN,
                        lambda line: line.lstrip().startswith("CTM bin 02"),
                        _bpch_to_cubes,
                        priority=7)

_BPCH_spec_le = format_picker.FormatSpecification(
                        'Binary Punch File (BPCH) v2 little-endian',
                        _LEADING_LINE_BIN_LE,
                        lambda line: line.lstrip().startswith("CTM bin 02"),
                        _bpch_to_cubes_le,
                        priority=6)

ifileformats.FORMAT_AGENT.add_spec(_BPCH_spec)
ifileformats.FORMAT_AGENT.add_spec(_BPCH_spec_le)

# wrappers for Iris's loading and saving cubes functions
load = iris.load
load_cube = iris.load_cube
load_cubes = iris.load_cubes
load_raw = iris.load_raw
#load_strict = iris.load_strict    # depreciated

save = iris.save

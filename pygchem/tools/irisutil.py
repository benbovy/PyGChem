# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2014 Benoît Bovy
# see license.txt for more details
#

"""
A bunch of utility functions based on the SciTools's Iris package.
"""

import numpy as np
import iris

from pygchem import grid

#-----------------------------------------------------------------------------
# Generic functions not specific to GEOS-Chem datasets
#-----------------------------------------------------------------------------


def set_dim_order(cube, coord_names):
    """
    Set the dimensions of `cube` in the order
    defined by `coord_names`, i.e., a list that must
    contain all of the cube's coordinate names.

    Change the order of the cube dimensions in-place.

    Warning: calling this method will trigger any
    deferred loading, causing the cube’s data array
    to be loaded into memory.
    """
    ordered_dims = [cube.coord_dims(cube.coord(cn))[0]
                    for cn in coord_names]
    cube.transpose(new_order=ordered_dims)


def get_dim_from_coord1d(cube, coord_name):
    """
    Get from `cube` the dimension of the
    1-dimensional coordinate specified by
    `coord_name` (string).

    Returns the dimension number of the coordinate.
    Returns None if the coordinate is not dimensional.

    Raises :err:`iris.exceptions.CoordinateMultiDimError`
    if the coordinate is not 1-dimensional.
    """
    c = cube.coord(coord_name)
    cdims = cube.coord_dims(c)

    if len(cdims) > 1:
        raise iris.exceptions.CoordinateMultiDimError(
            "'{}' coordinate must be 1-dimensional"
            .format(coord_name)
        )

    elif not len(cdims):
        return None

    return cdims[0]


def permute_dim_aux_coords(cube, dim_coord, aux_coord):
    """
    Permute a dimension coordinate with
    an auxilliary coordinate in a cube.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        the cube to process (changes will be made in-place).
    dim_coord : string or :class:`iris.coords.DimCoord`
        (name of) the dimension coordinate to permute.
    aux_coord : string or :class:`iris.coords.AuxCoord`
        (name of) the auxilliary coordinate to permute.
        must contain the same cube dimension than `dim_coord`
        and must satisfy the criteria required for dimensional
        coordinates (1D, numeric, and strictly monotonic).
    """
    if not isinstance(dim_coord, iris.coords.DimCoord):
        dim_coord = cube.coord(name=dim_coord)
    if not isinstance(dim_coord, iris.coords.AuxCoord):
        aux_coord = cube.coord(name=aux_coord)

    if cube.coord_dims(aux_coord) != cube.coord_dims(dim_coord):
        raise iris.exceptions.CoordinateNotRegularError(
            "Cannot permute '{}' and '{}' coordinates"
            .format(dim_coord.name(), aux_coord.name())
        )
    else:
        dim = cube.coord_dims(dim_coord)[0]

    new_dim_coord = iris.coords.DimCoord.from_coord(aux_coord)
    new_aux_coord = iris.coords.AuxCoord.from_coord(dim_coord)

    cube.remove_coord(dim_coord)
    cube.remove_coord(aux_coord)
    cube.add_dim_coord(new_dim_coord, dim)
    cube.add_aux_coord(new_aux_coord, dim)


def remove_dim_aux_coords(cube, dimensions):
    """
    Remove in `cube` all auxilliary coordinates
    that contain `dimensions`.
    """
    for coord in cube.coords(dimensions=dimensions):
        if isinstance(coord, iris.coords.AuxCoord):
            cube.remove_coord(coord)


#-----------------------------------------------------------------------------
# Misc. fix functions
#-----------------------------------------------------------------------------

def fix_cube_attributes_vals(cube):
    """
    Iris issue ? when saving a cube to a netCDF file,
    a TypeError is raised (by netCDF4) if there are
    attributes with bool values.

    Convert it to int.
    """
    for key, val in cube.attributes.items():
        if type(val) == bool:
            cube.attributes[key] = int(val)


def fix_nd49_pressure(pedges_cube):
    """
    GEOS-Chem Bug ?? PSURF_PEDGE-$ is missing one vertical level (e.g., using
    the GEOS5_47L model, it has 47 levels but should have 48 levels).

    Workaround: Add one level on top of the grid and assign low, arbitrary
    pressure values at that level.

    `pedges` must have 4-dimensions
    (time, longitude, latitude, model_level_number)

    Returns a new cube with one additional top level.

    Procedure:
    - extract the top level of the pressure cube (e.g., level 47 for GEOS5_47L)
    - duplicate the cube
    - replace the vertical coordinate value for the duplicated cube
      (e.g., assign 48 for GEOS5_47L)
    - assign new, arbitrary low pressures to the duplicated cube
      (e.g., press_at_level_48 = 0.1 * pressure_at_level_47)
    - re-assemble the cubes - split/slices and re-merge - to obtain a
      new cube with 48 levels
    """
    pedges_cube_slices = list(pedges_cube.slices(['time',
                                                  'longitude',
                                                  'latitude']))
    pedges_cube_top = pedges_cube_slices[-1]

    pedges_cube_overtop = pedges_cube_top.copy()

    l_coord = iris.coords.DimCoord(np.array([48]),
                                   standard_name='model_level_number')
    pedges_cube_overtop.replace_coord(l_coord)

    low_pressures = pedges_cube_top.data * 0.1
    pedges_cube_overtop.data = low_pressures

    pedges_cube_slices.append(pedges_cube_overtop)
    # works only with both concatenate and merge
    # (iris issue ? normally only merge is needed!)
    pedges_cube_new = iris.cube.CubeList(pedges_cube_slices)
    pedges_cube_new = pedges_cube_new.concatenate().merge()[0]
    pedges_cube_new.transpose(new_order=[1, 2, 3, 0])

    return pedges_cube_new


#-----------------------------------------------------------------------------
# Conversion functions
#-----------------------------------------------------------------------------

def ppbC_2_ppbv(cube):
    """
    Convert to ppbv units for hydrocarbon tracers that
    have ppbC units.

    ppbC = parts per billion carbon
         = ppbv * number of carbon atoms in the tracer molecule
    """
    is_ppbC = cube.attributes.get('no_udunits2') == 'ppbC'
    is_hydrocarbon = cube.attributes.get('hydrocarbon')

    if is_hydrocarbon and is_ppbC:
        carbon_weight = cube.attributes.get('carbon_weight')
        cube.data = cube.data / carbon_weight
        cube.units = 'ppbv'


#-----------------------------------------------------------------------------
# Functions for computing grid-related data and other additional data
#-----------------------------------------------------------------------------

def gcgrid_2_coords(model_name, model_resolution,
                    region_box=None):
    """
    Get the X,Y,Z coordinates of the GEOS-Chem grid given
    by `model_name` and `model_resolution`.

    Parameters
    ----------
    model_name : string
        name of a GEOS-Chem grid model supported by pygchem
    model_resolution : (float, float)
        horizontal grid resolution (lon, lat), in degrees
    region_box: (int, int, int, int, int, int) or None
        grid indices of the 3D region box of interest
        (imin, imax, jmin, jmax, lmin, lmax).
        i: longitude, j: latitude, l: vertical levels

    Returns
    -------
    i_coord, j_coord, l_coord : :class:`iris.coords.DimCoord`
        Iris dimensional coordinates objects
        (longitude, latitude, model_level_number)

    """
    g = grid.CTMGrid.from_model(model_name,
                                resolution=model_resolution)

    lon_points, lat_points = g.lonlat_centers

    elon, elat = g.lonlat_edges
    lon_bounds = np.column_stack((elon[0:-1], elon[1:]))
    lat_bounds = np.column_stack((elat[0:-1], elat[1:]))

    levels = np.arange(1, g.Nlayers + 1)

    if region_box is not None:
        imin, imax, jmin, jmax, lmin, lmax = region_box

        lon_points = lon_points[imin-1:imax]
        lon_bounds = lon_bounds[imin-1:imax]

        lat_points = lat_points[jmin-1:jmax]
        lat_bounds = lat_bounds[jmin-1:jmax]

        levels = np.arange(lmin, lmax + 1)

    spherical_geocs = iris.coord_systems.GeogCS(
        iris.analysis.cartography.DEFAULT_SPHERICAL_EARTH_RADIUS
    )

    i_coord = iris.coords.DimCoord(lon_points,
                       bounds=lon_bounds,
                       standard_name="longitude",
                       var_name="longitude",
                       units="degrees_east",
                       coord_system=spherical_geocs)

    j_coord = iris.coords.DimCoord(lat_points,
                       bounds=lat_bounds,
                       standard_name="latitude",
                       var_name="latitude",
                       units="degrees_north",
                       coord_system=spherical_geocs)

    l_coord = iris.coords.DimCoord(levels,
                       standard_name="model_level_number",
                       var_name="model_level_number",
                       units="1")

    return i_coord, j_coord, l_coord


def get_altitude_coord(gridbox_heights, topography):
    """
    Compute the altitude (elevation a.s.l) coordinate
    from grid box heights and topography.

    Parameters
    ----------
    gridbox_heights : :class:`iris.cube.Cube`
        a 1-d, 2-d, 3-d or 4-d cube.
        The cube must have at least 'longitude',
        'latitude' and 'model_level_number' coordinates.
        'longitude' and 'latitude' can be either
        1-dimensional or scalar.
        'model_level_number' must be 1-dimensional and
        should include all vertical levels.
    topography : :class:`iris.cube.Cube`
        a 2D cube, which must have 'longitude' and
        'latitude' dimension coordinates. `topography`
        must at least cover the horizontal extent of
        `gridbox_heights` and have a compatible horizontal
        grid (i.e., coordinates points and bounds must match).

    Returns
    -------
    a :class:`iris.AuxCoord` object
        An auxilliary coordinate that contains
        space and time-dependent altitude values for
        the points and bounds of each grid cell.
        The coordinate has the same dimensions and
        the same units than the `gridbox_heights` cube.
        The coordinate also contains attributes of
        the `topography` cube.

    Examples
    --------
    >>> ctm_cubes = iris.load('ctm.bpch')
    >>> box_heights = ctm_cubes.extract_strict('BXHEIGHT_BXHGHT-$')
    >>> nox_tracer = ctm_cubes.extract_strict('NOx_IJ-AVG-$')
    >>> global_topography = iris.load_cube('dem_GEOS57_2x2.5_awm.nc')
    >>> altitude_coord = get_altitude_coord(box_heights,
    ...                                     global_topography)
    >>> nox_tracer.add_aux_coord(altitude_coord,
                                 data_dims=range(0, nox_tracer.ndim))

    """

    # extract the region of topography that cooresponds to
    # the extent of gridbox_heights
    gbh_lat = gridbox_heights.coord('latitude').points
    gbh_lon = gridbox_heights.coord('longitude').points

    subset = iris.Constraint(latitude=gbh_lat,
                             longitude=gbh_lon)
    topography_subset = topography.extract(subset)

    # ensure same units for gridbox_heights and topography
    topography_subset.convert_units(gridbox_heights.units)

    # ensure proper array broadcasting when computing altitude values
    topo_lon_coord = topography_subset.coord('longitude')
    topo_lat_coord = topography_subset.coord('latitude')
    topo_dim_coords = topography_subset.dim_coords

    if topography_subset.ndim == 0:
        base_level_altitude = topography_subset.data

    else:
        gbh_lon_dim = get_dim_from_coord1d(gridbox_heights,
                                           'longitude')
        gbh_lat_dim = get_dim_from_coord1d(gridbox_heights,
                                           'latitude')

        if (topo_lon_coord in topo_dim_coords
            and topo_lat_coord in topo_dim_coords):

            set_dim_order(topography_subset, ['longitude', 'latitude'])
            dim_map = (gbh_lon_dim, gbh_lat_dim)

        elif topo_lon_coord in topo_dim_coords:
            dim_map = (gbh_lon_dim,)
        else:
            dim_map = (gbh_lat_dim,)

        base_level_altitude = iris.util.broadcast_to_shape(
            topography_subset.data.copy(),
            gridbox_heights.shape,
            dim_map)

    level_dim = get_dim_from_coord1d(gridbox_heights,
                                     'model_level_number')

    # calculate altitude points and bounds values
    altitude_ubnd = np.cumsum(gridbox_heights.data, axis=level_dim) + \
                    base_level_altitude
    altitude_lbnd = altitude_ubnd - gridbox_heights.data
    altitude_lbnd[..., 1:] = altitude_ubnd[..., :-1]  # force contiguous bounds
    altitude_points = (altitude_lbnd + altitude_ubnd) / 2.

    altitude_bounds = np.concatenate((altitude_lbnd[...,np.newaxis],
                                      altitude_ubnd[...,np.newaxis]),
                                     axis=-1)

    # create and return the iris coordinate object
    altitude_coord = iris.coords.AuxCoord(
        altitude_points,
        bounds=altitude_bounds,
        standard_name='altitude',
        units=gridbox_heights.units,
        attributes=topography.attributes.copy()
    )

    return altitude_coord


def compute_cell_volumes(cube, lon_coord='longitude', lat_coord='latitude',
                         height_coord='altitude', as_new_cube=True):
    """
    Compute the volume of each cell of a grid
    defined by its (spherical) coordinates.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube` object
        a cube that contain the latitude, longitude and
        height coordinates.
    lon_coord : string or :class:`iris.coords.Coord` object
        longitude horizontal coordinate (1-dimensional
        or scalar), expressed either in degrees or radians
        (must have bounds).
    lat_coord : string or :class:`iris.coords.Coord` object
        latitude horizontal coordinate (1-dimensional
        or scalar), expressed either in degrees or radians
        (must have bounds).
    height_coord : string or :class:`iris.coords.Coord` object
        altitude/height vertical coordinate.
        This coordinate must also have bounds.
        It can be either 1-dimensional or multi-dimensional.
    as_new_cube : bool
        if True, returns a new cube'. Otherwise it adds
        to `cube` the computed volumes values as a new
        auxilliary coordinate.

    Raises
    ------
    :err:`iris.exceptions.CoordinateMultiDimError`
        if latitude and/or longitude coordinates
        are multi-dimensional.
    :err:`iris.exceptions.CoordinateNotRegularError`
        if latitude and/or longitude coordinates
        don't use a spherical coordinate system.

    Notes
    -----
    the altitude coordinate should be orthogonal to lat/lon
    coordinates.

    Volumes are calculated for each cell as:

    .. math::

        \frac{1}{3} \left(\theta_{0} - \theta_{1}\right)
        \left(\rho_{0}^{3} \cos{\left (\phi_{0} \right )} -
        \rho_{0}^{3} \cos{\left (\phi_{1} \right )} -
        \rho_{1}^{3} \cos{\left (\phi_{0} \right )} +
        \rho_{1}^{3} \cos{\left (\phi_{1} \right )}\right)

    where :math:`\rho_0` and :math:`\rho_1` are the altitude
    bounds (+ earth radius), :math:`\phi_0` and :math:`\phi_1`
    are the latitude bounds and :math:`\theta_0` and :math:`\theta_1`
    are the longitude bounds.

    Uses earth radius from the lat/lon coordinates,
    if present and spherical. Defaults to
    iris.analysis.cartography.DEFAULT_SPHERICAL_EARTH_RADIUS.

    """
    # get coordinate objects
    if not isinstance(lon_coord, iris.coords.Coord):
        lon_coord = cube.coord(lon_coord)
    if not isinstance(lat_coord, iris.coords.Coord):
        lat_coord = cube.coord(lat_coord)
    if not isinstance(height_coord, iris.coords.Coord):
        height_coord = cube.coord(height_coord)

    # verify that lon_coord and lat_coord are 1-dimensional
    if lon_coord.ndim > 1 or lat_coord.ndim > 1:
        raise iris.exceptions.CoordinateMultiDimError(
            "multi-dimensional longitude and/or latitude "
            "coordinates are not supported"
        )

    # copy coordinates and convert units
    lon_coord_rad = lon_coord.copy()
    lon_coord_rad.convert_units('radians')
    lat_coord_rad = lat_coord.copy()
    lat_coord_rad.convert_units('radians')
    height_coord_m = height_coord.copy()
    height_coord_m.convert_units('m')

    out_units = height_coord.units**3

    # altitude -> radius (in meters)
    if lon_coord.coord_system is not None:
        err_msg = ("latitude and longitude coordinates don't use "
                   "a spherical coordinate system")
        cs = lon_coord.coord_system
        if not isinstance(cs, iris.coord_systems.GeogCS):
            raise iris.exceptions.CoordinateNotRegularError(err_msg)
        if cs.semi_major_axis != cs.semi_minor_axis:
            raise iris.exceptions.CoordinateNotRegularError(err_msg)
        earth_radius = cs.semi_major_axis

    if lat_coord.coord_system != lon_coord.coord_system:
        raise iris.exceptions.CoordinateNotRegularError(
            "latitude and longitude coordinates don't use "
            "the same coordinate system"
        )

    else:
        earth_radius = iris.analysis.cartography.DEFAULT_SPHERICAL_EARTH_RADIUS

    height_bounds = height_coord_m.bounds + earth_radius

    # broadcast lat_coord and/or lon_coord bounds if needed
    if lon_coord.points.size > 1:
        lon_dim = cube.coord_dims(lon_coord)[0]
        dim_map = (lon_dim, -1)
        lon_bounds = iris.util.broadcast_to_shape(lon_coord_rad.bounds,
                                                  height_bounds.shape,
                                                  dim_map)
    else:
        lon_bounds = lon_coord_rad.bounds

    if lat_coord.points.size > 1:
        lat_dim = cube.coord_dims(lat_coord)[0]
        dim_map = (lat_dim, -1)
        lat_bounds = iris.util.broadcast_to_shape(lat_coord_rad.bounds,
                                                  height_bounds.shape,
                                                  dim_map)
    else:
        lat_bounds = lat_coord_rad.bounds

    # calculate grid-cell volumes
    phi_0, phi_1 = lat_bounds[..., 0], lat_bounds[..., 1]
    theta_0, theta_1 = lon_bounds[..., 0], lon_bounds[..., 1]
    rho_0, rho_1 = height_bounds[..., 0], height_bounds[..., 1]

    volume = rho_0**3 * np.cos(phi_0) - \
             rho_0**3 * np.cos(phi_1) - \
             rho_1**3 * np.cos(phi_0) + \
             rho_1**3 * np.cos(phi_1)
    volume *= 1./3. * (theta_0 - theta_1)

    # convert volume units
    m = iris.unit.Unit('m^3')
    volume = m.convert(volume, out_units)

    # either return a new iris cube or add aux coordinate
    v_long_name = 'grid-cell_true_volume'

    if as_new_cube:
        vol_cube = cube.copy()
        vol_cube.data = volume
        vol_cube.units = out_units
        vol_cube.long_name = v_long_name
        vol_cube.standard_name = None
        return vol_cube

    else:
        cube.add_aux.coord(iris.coords.AuxCoord(
            volume,
            long_name=v_long_name,
            units=out_units,
            data_dims=range(0, height_coord.ndim)))


def compute_tracer_columns(mixing_ratio, n_air, z_coord,
                           cell_heights=None, weights=None,
                           units='count/cm2'):
    """
    Compute tracer total columns.

    Parameters
    ----------
    mixing_ratio : :class:`iris.cube.Cube` object
        Tracer mixing ratio.
    n_air : :class:`iris.cube.Cube` object
        Number density of air.
    z_coord : string or :class:`iris.coords.Coord` object
        (Name of) the vertical levels coordinate.
    cell_heights : :class:`iris.cube.Cube` object or None
        Grid cell heights. Must be provided if `z_coord`
        can't be used to compute the height of each cell.
    weights : :class:`iris.cube.Cube` object
        vertical cell weights.
    units : string or :class:iris.unit.`Unit` object
        Units in which the results are returned.

    Returns
    -------
    A :class:`iris.cube.Cube` object with the computed
    total columns values.

    Notes
    -----
    - This function computes total columns.
      To compute partial columns, define `weights`
      or reduce the cubes first.
    - Partial columns are first calculated for each grid cell
      as follows:
        pcol = `mixing_ratio` * `N_air` * `cell_heights`
      total columns are the sum of pcol over `z_coord`.
    - Aware of units (e.g., ppbv, ppmv for `mixing_ratio`...).
      specified in the input cubes metadata.
    - For output units, 'count' is used instead of 'molec' (the latter
      unit is not supported by udunits2), but the values are unchanged.

    """
    if not isinstance(z_coord, iris.coords.Coord):
        z_coord = mixing_ratio.coord(name=z_coord)

    if weights is None:
        weights = 1.

    if cell_heights is None:
        if not z_coord.has_bounds():
            raise iris.exceptions.CoordinateNotRegularError(
                'z_coord {} must have bounds'.format(z_coord.name)
            )
        if z_coord.units.symbol not in ['m', 'km']:
            raise iris.exceptions.CoordinateNotRegularError(
                "'m' or 'km' units required for z_coord, "
                "found {} ({})".format(z_coord.units.symbol, z_coord.name)
            )

        map_height_dim = mixing_ratio.coord_dims(z_coord)
        heights_data = z_coord.bounds[:,1] - z_coord.bounds[:,0]
        heights_data = iris.util.broadcast_to_shape(heights_data,
                                                    mixing_ratio.shape,
                                                    map_height_dim)
        cell_heights = mixing_ratio.copy()
        cell_heights.data = heights_data
        cell_heights.units = z_coord.units

    # avoid missing values (convert to zero)
    mixing_ratio = mixing_ratio.copy()
    n_air = n_air.copy()
    mixing_ratio.data = np.nan_to_num(mixing_ratio.data)
    n_air.data = np.nan_to_num(n_air.data)

    gridbox_partial_columns = mixing_ratio * n_air
    gridbox_partial_columns *= cell_heights * weights

    total_columns = gridbox_partial_columns.collapsed(z_coord,
                                                      iris.analysis.SUM)
    total_columns.convert_units(units)
    total_columns.long_name = "_".join([mixing_ratio.name(), 'columns'])
    total_columns.attributes.update(mixing_ratio.attributes)

    return total_columns


#-----------------------------------------------------------------------------
# Regridding functions
#-----------------------------------------------------------------------------

def _compute_exchange_vertical(src_z_coord, dst_z_coord):
    """
    Compute the exchange vertical grid needed for
    vertical regridding.

    The cells of the exchange grid are defined by the
    intersections of the cells of the source grid and the
    the cells of the destination grid.

    Returns a :class:`iris.coords.DimCoord` object.

    """
    temp = np.concatenate((src_z_coord.bounds,
                           dst_z_coord.bounds))
    temp = np.sort(temp, kind='mergesort', axis=None)
    temp = np.unique(temp)

    exchange_z_bounds = np.column_stack((temp[:-1], temp[1:]))
    exchange_z_points = (exchange_z_bounds[:, 0] + exchange_z_bounds[:, 1]) / 2.

    exchange_z_coord = iris.coords.DimCoord(
        exchange_z_points,
        bounds=exchange_z_bounds,
        standard_name=src_z_coord.standard_name,
        units=src_z_coord.units
    )

    return exchange_z_coord


def _map_exchange_vertical(src_or_dst_data,
                           src_or_dst_z_coord,
                           exchange_z_coord):
    """
    Map on the exchange vertical grid data values given
    on the source (or the destination) vertical grid.

    Returns only the data (numpy array).

    If a cell of the exchange grid doesn't overlap any
    of the source (destination grid), value of that cell
    is set to `np.nan`.

    """
    points = exchange_z_coord.points
    bounds = src_or_dst_z_coord.bounds

    # 2-d boolean array (rows are exchange cells
    # and cols are source/destination grid cells)
    # True if exchange cell match with src or dst cell
    is_in_cell = ((points[:,np.newaxis] > bounds[:,0]) &
                  (points[:,np.newaxis] < bounds[:,1]))

    # a cell of the exchange should have only one
    # corresponding cell in the source/destination grid
    if is_in_cell.sum(axis=1).max() > 1:
        raise ValueError("source/destination grid is not "
                         "compatible with the exchange grid")

    # indices map between exchange grid (1st array)
    # and source or destination grid (2st array)
    map_indices = np.nonzero(is_in_cell)

    # assign data
    map_data = np.empty_like(points) * np.nan
    map_data[map_indices[0]] = src_or_dst_data[map_indices[1]]

    return map_data


def regrid_conservative_vertical(src_cube, dst_grid_cube,
                                 src_data_type='intensive',
                                 z_coord_name='altitude',
                                 overlap_tol=0.95):
    """
    Conservative regridding of a vertical profile.

    Parameters
    ----------
    src_cube : :class:`iris.cube.Cube`
        the original profile to be regridded
    dst_grid_cube : :class:`iris.cube.Cube`
        defines the destination grid
    src_data_type : ('intensive', 'extensive')
        specifies whether the data of `src_cube`
        is an intensive field (i.e., quantities
        whose value changes with the grid cell size)
        or an extensive field (i.e., grid-size
        independent quantities).
    z_coord_name : string
        name of the vertical coordinate to use
        for regridding (must be present in both
        `src_cube` and `dst_grid_cube`).
    overlap_top : float
        cell overlap tolerance, i.e., a threshold
        of - relative [0, 1] - overlapping height
        between a cell of the destination grid
        and cells of the source grid, under which
        the value of the destination cell is set
        to undefined.

    Returns
    -------
    A new :class:`iris.cube.Cube` object
        a copy of `dst_grid_cube` with
        the regridded data of `src_cube`.
        Auxilliary / scalar coordinates and
        attributes of `src_cube` will be copied.

    Notes
    -----
    Both `src_cube` and `dst_grid_cube` must be
    1-dimensional and must have compatible
    vertical coordinates, i.e., 1-dimensional,
    with the same name `z_coord_name`, the same
    coordinate system, the same units
    and contiguous bounds.


    """
    # get or compute vertical coordinates of the source, destination
    # and exchange grids
    # TODO: check if vertical coordinates of src_cube and dst_cube
    # exist, are 1-dimensional and are compatibles (units,
    # coord system, contiguous bounds)
    src_z_coord = src_cube.coord(name=z_coord_name)
    dst_z_coord = dst_grid_cube.coord(name=z_coord_name)

    exchange_z_coord = _compute_exchange_vertical(src_z_coord,
                                                  dst_z_coord)

    # compute cell heights of each grid
    src_heights = src_z_coord.bounds[:,1] - src_z_coord.bounds[:,0]
    dst_heights = dst_z_coord.bounds[:,1] - dst_z_coord.bounds[:,0]
    exchange_heights = exchange_z_coord.bounds[:,1] - \
                       exchange_z_coord.bounds[:,0]

    # generate arbitray levels for the destination grid
    # (used for further aggregation - summing - of the exchange grid
    # cells onto the destination grid
    dst_levels = np.arange(1, dst_z_coord.points.size + 1)

    # map source data values, source heights and destination heights
    # (and destination levels) on the exchange grid
    src_heights_mapped = _map_exchange_vertical(src_heights,
                                                src_z_coord,
                                                exchange_z_coord)
    dst_heights_mapped = _map_exchange_vertical(dst_heights,
                                                dst_z_coord,
                                                exchange_z_coord)
    dst_levels_mapped = _map_exchange_vertical(dst_levels,
                                               dst_z_coord,
                                               exchange_z_coord)
    src_data_mapped = _map_exchange_vertical(src_cube.data,
                                             src_z_coord,
                                             exchange_z_coord)

    # compute interpolation weights and data values
    # on the exchange grid
    if src_data_type == 'intensive':
        exchange_iweights = exchange_heights / dst_heights_mapped
    elif src_data_type == 'extensive':
        exchange_iweights = exchange_heights / src_heights_mapped
    else:
        raise ValueError("invalid source data type: {}"
                         .format(src_data_type))

    exchange_data = exchange_iweights * src_data_mapped

    # finally get data values on the destination grid
    # (sum data of overlapping cells of the exchange grid
    # and apply the overlap tolerance)
    levels_mask = dst_levels[:,np.newaxis] == dst_levels_mapped[np.newaxis,:]

    dst_data = np.nansum((exchange_data[np.newaxis,:] * levels_mask),
                         axis=1)

    valid_heights = np.where(np.isnan(exchange_data), np.nan, exchange_heights)
    overlap_heights = np.nansum((valid_heights[np.newaxis,:] * levels_mask),
                                axis=1)
    overlap_heights_relative = overlap_heights / dst_heights
    nan_indices = np.nonzero(overlap_heights_relative < overlap_tol)
    dst_data[nan_indices] = np.nan

    # construct the cube of the regridded field
    # copy all scalar coordinates from the source cube
    dst_cube = dst_grid_cube.copy()
    dst_cube.data = np.array(dst_data)
    dst_cube.var_name = "_".join([src_cube.name(), "regridded"])
    dst_cube.long_name = "_".join([src_cube.name(), "regridded"])
    dst_cube.attributes.update(src_cube.attributes)
    dst_cube.attributes.update({'src_cube': src_cube.name(),
                                'grid_cube': dst_cube.name(),
                                'regrid_method': 'first-order conservative'})
    dst_cube.units = src_cube.units

    for c in src_cube.coords():
        if c.points.size == 1 and not dst_cube.coords(coord=c):
            dst_cube.add_aux_coord(c)

    return dst_cube

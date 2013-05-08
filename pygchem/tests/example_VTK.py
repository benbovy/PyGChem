# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 10:18:46 2012

@author: bovy
"""

import os, sys

import numpy as np
from tvtk.api import tvtk

from pygchem import diagnostics
from pygchem import grid as gchemgrid


ws = "/Volumes/data/03_geo/geoschem/GEOS-Chem-rundirs/4x5/geos5/standard"
bpch_fname = "ctm.bpch"

# open ctm bunch file
ctm_f = diagnostics.CTMFile.fromfile(bpch_fname)



# vtk spherical mesh

def latlonr2xyz(points):
    """
    convert earth spherical point coordinates (lat/lon/radius) to their
    cartesian equivalent (xyz)

    (x,y,z) cartesian coordinates are based on the following axes:
        - x-axis: center of the earth pointing towards 0 lat/long
        - y-axis: center of the earth pointing towards 0 lat / 90 lon
        - z-axis: center of the earth pointing towards north pole

    points: 2D array (Nx3), one row = one point, columns = (lon, lat, radius),
            lon/lat must be given in degrees, and radius in the unit of the
            catersian system

    return Nx3 2D array
    """

    points[:,0:2] *= pi/180.0
    points_xyz = empty_like(points)
    points_xyz[:,0] = points[:,2] * np.cos(points[:,0]) * np.cos(points[:,1])
    points_xyz[:,1] = points[:,2] * np.sin(points[:,0]) * np.cos(points[:,1])
    points_xyz[:,2] = points[:,2] * np.sin(points[:,1])


    return points_xyz


def xyz2latlonr(points):
    """
    convert earth spherical cartesian point coordinates (xyz) to their
    (lat/lon/radius) equivalent

    (x,y,z) cartesian coordinates are based on the following axes:
        - x-axis: center of the earth pointing towards 0 degree lat/long
        - y-axis: center of the earth pointing towards 0 degree lat / 90 lon
        - z-axis: center of the earth pointing towards north pole

    points: Nx3 2D array, one row = one point, columns = (x,y,z) coordinates

    return 2D array (Nx3), columns = (lon, lat, radius),
            lon/lat must are given in degrees, and radius in the unit of the
            catersian system
    """

    points_llr = empty_like(points)
    points_llr[:,2] = np.sqrt(np.sum(np.power(points, 2.0), axis=1))
    points_llr[:,1] = np.arctan2(points[:,2],
                        np.sqrt(np.sum(np.power(points[:,0:2], 2.0), axis=1)))
    points_llr[:,0] = np.arctan2(points[:,1],points[:,0])
    points_llr[:,0:2] *= 180.0/pi

    return points_llr


def smesh_points(lon,lat,lev,lev_scale,lev_offset):
    """
    return a list of N (x,y,z) points coordinates (as a Nx3 array)
    from longitude-latitude-level (or radius) N-sized 1D arrays,
    which defines a 3D spherical mesh
    Useful as input of a vtk structured grid

    order of points append to the output array:
        parallels > spheres layers (lon>lat>levels)

    lev_scale: scaling factor for level values (in radius units)
    lev_offset: (radius) value of 1st level, a value > 0 means that
                the mesh doesn't cover the center of the sphere

    call routines:
        latlong2xyz
    """

    # create array of point coords in lon-lat-lev, following specific order
    # each row is a grid point: dimension (N,3)
    points = np.column_stack((np.tile(np.tile(lon,lat.size),lev.size),
                              np.tile(np.repeat(lat,lon.size), lev.size),
                              np.repeat(lev*lev_scale+lev_offset, lon.size*lat.size))
                              )

    n_dims = (lon.size, lat.size, lev.size)

    # transform spherical coordinates (degrees for lat/lon) into cartesian coords
    points_xyz = latlonr2xyz(points)


    return n_dims, points_xyz


def smesh_cells(lon,lat,lev,lev_scale,lev_offset,lon_connect):
    """
    call routines:
        smesh_points
    """

    n_dims, points_xyz = smesh_points(lon,lat,lev,lev_scale,lev_offset)
    n_points = points_xyz.shape[0]

    # vertices indices for the 1st cell (warning, vertices order -> draw)
    ind_offset_list = [0, 1,
                       n_dims[0]+1, n_dims[0],
                       n_dims[0]*n_dims[1], n_dims[0]*n_dims[1]+1,
                       n_dims[0]*(n_dims[1]+1)+1, n_dims[0]*(n_dims[1]+1)
                      ]

    # build indices increments
    n_cells = n_points - ind_offset_list[-2]
    vertices = empty((8,n_cells),dtype='int')
    for i in xrange(0,8):
        vertices[i] = np.arange(ind_offset_list[i],
                                ind_offset_list[i]+n_cells,
                                1)
    # cells loop in longitude: indices correction or remove vertices
    loop_ind = arange(n_dims[0]-1,n_cells,n_dims[0])
    if(lon_connect == False):
        for i in [1,2,5,6]:
            vertices[i][loop_ind] -= n_dims[0]
    elif(lon_connect == True):
        vertices = np.delete(vertices, loop_ind, axis=1)

    # cells no-loop in latitude: remove vertices
    for i in xrange((n_dims[0]-1)*(n_dims[1]-1),n_cells,(n_dims[0]-1)*(n_dims[1]-1)):
        loop_ind = arange(i,i+n_dims[0]-1,1)
        vertices = np.delete(vertices, loop_ind, axis=1)


    return n_dims, points_xyz, vertices.transpose()


def diagnostics2vtk(diagnostics, mesh_type='spherical', geom='cell',
                    create_earth=True, exageration=2):
    """
    Export a list of diagnostics to VTK.
    
    Agruments
    ---------
    diagnostics : iterable
        sequence of CTM diagnostics
    mesh_type : string
        `spherical` or `projected`
    create_earth : bool
        generate VTK files for continents (and earth if mesh_type is spherical)
    exageration : float
        scale factor applied to vertical levels for easy visualization
    """
   
    # get scale factor for 
    # TODO : get vertical grid type (GEOS5, GEOS5 reduced...) from diagnostics
    atm_thickness = gchemgrid.c_km_geos5_r.max() * 1e3 
    



# vtk atmosphere thickness vs. earth radius: height scale calculation

ratio_earth_atm = 2.0

earth_radius_scaled = ratio_earth_atm * atm_thickness


# vtk file basemap (sphere, continents)
globe_src = tvtk.SphereSource(radius=earth_radius_scaled-1e-3*earth_radius_scaled,
                              lat_long_tessellation=True,
                              phi_resolution=gchemgrid.c_lat_4x5.size,
                              theta_resolution=gchemgrid.c_lon_4x5.size)
continents_src = tvtk.EarthSource(on_ratio=1, radius=earth_radius_scaled)

writer = tvtk.XMLPolyDataWriter(input=globe_src.output,
                        file_name=os.path.join(run_dir,"vtk","globe.vtp"))
writer.write()
writer = tvtk.XMLPolyDataWriter(input=continents_src.output,
                        file_name=os.path.join(run_dir,"vtk","continents.vtp"))
writer.write()


# vtk meshes
tfield = filter_results[0].values

lev_offset = earth_radius_scaled + 1e-3 * earth_radius_scaled
lev_scale = 1.0  #atm_thickness / tfield.shape[2]


# cell-centered point-based mesh
sdims, spts = smesh_points(gchemgrid.c_lon_4x5,
                     gchemgrid.c_lat_4x5,
                     gchemgrid.c_km_geos5_r * 1e3, #arange(tfield.shape[2]),
                     lev_scale,lev_offset)

sgrid = tvtk.StructuredGrid(dimensions=sdims)
sgrid.points = spts

sgrid.point_data.scalars = np.ravel(tfield.transpose())
sgrid.point_data.scalars.name = filter_results[0].full_name

writer = tvtk.XMLStructuredGridWriter(input=sgrid,
                        file_name=os.path.join(run_dir,"vtk","test.vts"))
writer.write()


# cell-based mesh (use unstructured mesh, polydata is for surfaces and not
# volumes)
sdims, spts, vertices = smesh_cells(gchemgrid.e_lon_4x5,
                          gchemgrid.e_lat_4x5,
                          gchemgrid.e_km_geos5_r * 1e3,
                          lev_scale,lev_offset,True)

sgrid = tvtk.UnstructuredGrid(points=spts)
tet_type = tvtk.Hexahedron().cell_type # cells = hexahedrons
sgrid.set_cells(tet_type, vertices)


ijavg = ctm_f.filter(category=field_category)
for tracer in ijavg:
    print "export diagnostic %s to VTK" % tracer.full_name
    vals_array = tracer.values * units_conv[mr_units]
    #ar = sgrid.cell_data.scalars = np.ravel(vals_array.transpose())
    #sgrid.cell_data.scalars.name = tracer.full_name
    ar = sgrid.cell_data.add_array(np.ravel(vals_array.transpose()))
    sgrid.cell_data.get_array(ar).name = tracer.full_name

outfile = os.path.join(run_dir,"vtk","test_ijavg_cell.vtu")
if os.path.exists(outfile):
    os.remove(outfile)
writer = tvtk.XMLUnstructuredGridWriter(input=sgrid,
                        file_name=outfile)
writer.write()

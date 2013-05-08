# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 10:18:46 2012

@author: bovy

Routines for plotting CTM diagnostics on 2D maps.
This is a preliminar code for features that will be further part of
the pygchem library.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from pygchem import grid as gchemgrid


def map_level_avg(diagnostic, bottom_level, top_level, **kwargs):
    """
    Create a 2D map of a 3D diagnostic by averaging between `bottom_level` and
    `top_level` (top_level not included).
    """
    
    # res2grids = ... make dict or list of correspondence instead of if blocks
    if diagnostic.resolution == (5, 4):
        lon = gchemgrid.c_lon_4x5
        lat = gchemgrid.c_lat_4x5
    
    scalar = diagnostic.values[:,:,bottom_level:top_level].mean(axis=2)
    scalar = scalar.transpose() * diagnostic.scale
    
    title = "%s for %s (Avg from L%d-%d)" % (diagnostic.full_name,
                                             diagnostic.times[0],
                                             bottom_level,
                                             top_level)
    
    return scalar_map2D(lat, lon, scalar, units=diagnostic.unit, title=title,
                        **kwargs)
    

def scalar_map2D(lat, lon, scalar, units="", title="", projection="cyl",
                 clb_loc="bottom", ):
    """
    Create a 2D map using basemap.
    
    Arguments
    ---------
    clb_loc : location of the colorbar ("top", "bottom", "left", "right")
    
    Return
    ------
    map2D : object
        The Basemap object
    plot : object
        The pcolormesh object
    clb : object
        The colorbar object
    """
    
    bmap = Basemap(projection=projection, lon_0=0.0, resolution='c')
    
    x, y = bmap(*np.meshgrid(lon, lat))
    plot = bmap.pcolormesh(x, y, scalar, cmap=plt.cm.Blues)
    
    # draw coastlines.
    bmap.drawcoastlines(color=(0.2,0.2,0.2), linewidth=0.6)
    # draw a line around the map region.
    bmap.drawmapboundary()
    # draw parallels and meridians.
    bmap.drawparallels(np.arange(-60.,90.,30.), labels=[1,0,0,0],
                        color=(0.4,0.4,0.4), linewidth=(0.6))
    bmap.drawmeridians(np.arange(-180.,180.,60.), labels=[0,0,0,1],
                        color=(0.4,0.4,0.4), linewidth=(0.6))
    # draw colorbar
    clb = bmap.colorbar(mappable=plot, location=clb_loc,
                         size='5%', pad='15%')
    clb.ax.set_ylabel("[%s]" % units, rotation='horizontal',
                      labelpad=10.0)
    clb.ax.yaxis.set_label_position('right')
    
    plt.title(title)
    plt.show()
    
    return bmap, plot, clb
    
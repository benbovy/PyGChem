# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 10:18:46 2012

@author: bovy

A simple example for mapping CTM diagnostics using matplotlib/basemap
"""

import os
import datetime

import temp.diagnostics as gdiag
import temp_map2D


ws = "/Volumes/data/03_geo/geoschem/GEOS-Chem-rundirs/4x5/geos5/standard"
bpch_fname = "ctm.bpch"

# open ctm bunch file
ctm_f = gdiag.CTMFile.fromfile(os.path.join(ws, bpch_fname))

# select diagnostic(s)
diagnostics = ctm_f.filter(name="Ox", category="IJ-AVG-$",
                           time=datetime.datetime(2005,07,1))

# create 2D map(s) by averaging over a level interval
for diag in diagnostics:
    temp_map2D.map_level_avg(diag, 0, 10)

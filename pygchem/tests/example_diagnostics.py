# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 14:50:17 2013

@author: bovy

Examples using the diagnostics module of pygchem.
"""

import os
import datetime

import numpy as np

import pygchem.diagnostics as gdiag


ws = "/Volumes/data/03_geo/geoschem/GEOS-Chem-rundirs/4x5/geos5/standard"
bpch_fname = "ctm.bpch"

# import a CTM file (only binary punch file format is supported,
# netCDF will follow)
ctm_f = gdiag.CTMFile.fromfile(os.path.join(ws, bpch_fname))

# diagnostics metadata (contents of diaginfo.dat and tracerinfo.dat) are
# automatically imported. Access to metadata via an instance of the
# 'Diagnostics' class connected to the CTMFile instance
ctm_f.diagnostics             # returns a 'Diagnostics' instance
ctm_f.diagnostics.categories  # returns a dict with categories (diaginfo)
ctm_f.diagnostics.diagnostics # returns a dict with diagnostics (tracerinfo)

# Possible to add/remove/modify diagnotics metadata at runtime...
new_diags = gdiag.Diagnostics()
new_diags.add_diagnostic_field('newall', 'newall_val', True)

# access to datablocks in the CTM file is provided by : 
ctm_f.datablocks              # returns a list of datablock (as 'Datablock'
                              # instances)

# datablock simple filtering (by name, number, category and/or times) :
ijavg = ctm_f.filter(category="IJ-AVG-$")   # select all IJ-AVG-$ datablocks
d_050701 = ctm_f.filter(time=datetime.datetime(2005,07,1))  # time selection
Ox_avg_050701 = ctm_f.filter(name="Ox",                     # Ox tracer...
                             category="IJ-AVG-$",
                             time=datetime.datetime(2005,07,1))[0]  

# datablock advanced filtering (TODO)

# datablock header (examples) and values
Ox_avg_050701.index
Ox_avg_050701.number
Ox_avg_050701.name
Ox_avg_050701.full_name
Ox_avg_050701.molecular_weight
Ox_avg_050701.unit

Ox_avg_050701.values      # returns a numpy array


# create new datablocks
default_diags = gdiag.Diagnostics()    # create a new Diagnotics instance
                                       # (here based on 'default' tracerinfo
                                       # and diaginfo)

d_start = datetime.datetime(2000,1,1)  # start and end times assigned to 
d_end = datetime.datetime(2002,1,1)    # datablock

new_db = gdiag.DataBlock(1, 'EW-FLX-$', (d_start, d_end), 
                         diagnostics=default_diags,
                         values=np.zeros_like(Ox_avg_050701.values))

# append the new datablock to a new 'CTMFile' instance and write to a bunch
# binary file (useful for creating restart files)
new_ctm_f = gdiag.CTMFile(diagnostics=default_diags)
new_ctm_f.append_datablock(new_db)
new_ctm_f.save(os.path.join(ws, "new_ctm.bpch"), overwrite=True)

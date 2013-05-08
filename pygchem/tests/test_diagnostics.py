# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 14:50:17 2013

@author: bovy
"""

import os
import sys
#sys.path.append("/home/bovy/code_projects/python_geoschem/PyGChem")
import pygchem.diagnostics as gdiag
import datetime

#import cProfile

dstart = datetime.datetime(2000,1,1)
dend = datetime.datetime(2002,1,1)

standard_diags = gdiag.Diagnostics()

standard_diags.add_diagnostic_field('newall', 'newall_val', True)
standard_diags.add_diagnostic_field('new5', 'new5_val', False, 201)

db = gdiag.DataBlock(1, 'EW-FLX-$', (dstart, dend),
                     diagnostics=standard_diags)

db2 = gdiag.DataBlock(2, 'EW-FLX-$', (dstart, dend), name='testtracer')


ws = "/Volumes/data/03_geo/geoschem/GEOS-Chem-rundirs/4x5/geos5/standard"
bpch_fname = "ctm.bpch"

ctm = gdiag.CTMFile.fromfile(os.path.join(ws, bpch_fname))

#ctm.save(os.path.join(ws, "test.bpch"), overwrite=True)

#ctm2 = gdiag.CTMFile.fromfile(os.path.join(ws, "test.bpch"))
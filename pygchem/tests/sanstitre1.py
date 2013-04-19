# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 18:01:16 2012

@author: bovy
"""

import os
import sys
sys.path.append("/home/benbovy/code_projects/python_geoschem/PyGChem")
import pygchem.globchem as gl
import pygchem.io.globchem as iogl

import cProfile


ws = "/home/data/03_geo/geoschem/GEOS-Chem-rundirs/4x5/geos5/standard/"
fn_globchem = os.path.join(ws, "globchem.dat")


mb_groups, species, reac_kn, reac_ph = iogl.read_globchem_dat(fn_globchem)

moh = gl.Species.from_dict(species['MOH'])

newchem = {'MOH': moh}
# python set built-in type (from python 2.4): for 
compared =  set(newchem['MOH'].to_dict()) ^ set(newchem['MOH'].__dict__)

cProfile.run('iogl.read_globchem_dat(fn_globchem)')

#gfile = gl.Globchem()
#gfile.from_file(fn_globchem)




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

# test io.globchem.read_globchem_dat
(header, mb_groups,
 species, reac_kn, reac_ph) = iogl.read_globchem_dat(fn_globchem)
#cProfile.run('iogl.read_globchem_dat(fn_globchem)')  # benchmark

# test create a species
moh = gl.Species.from_dict(species['MOH'])
newchem = {'MOH': moh}
compared =  set(newchem['MOH'].to_dict()) ^ set(newchem['MOH'].__dict__)
     # python set built-in type (from python 2.4): for

# test globchem
gc = gl.Globchem()
for s in species.values():
    gc.add_species(gl.Species.from_dict(s))

del species['CO2']['name']
co2 = gl.Species.from_dict(species['CO2'])
# give a UserWarning: species attribute 'name' not in the dictionary...

gc.remove_species(['MOH', 'GLYC'], raise_error=False)

gc.copy_species('H2O', ['H2O_A', 'H2O_B'], keep_original=False,
                name=["water vapor A", "water vapor B"],
                formula=["H2O A", "H2O B"])
print gc.species['H2O_A'].to_dict()


gc2 = gl.Globchem.from_smvgear2(fn_globchem)

print gc2.filter_mb_groups_bystatus('A')
# returns: ['CAR', 'SUL', 'CHL', 'BRO', 'NH4', 'FLO', 'NO3']
gc2.add_mb_group('ADD', value=1)
gc2.remove_mb_group('HYD')

test_expr = lambda formula, status: formula in ("H2O", "CO2") and status == "D"
print gc2.filter_species(test_expr)
# returns: ['CO2']

gc2.remove_species("HO2")
# give a UserWarning: kinectic reaction(s) involving species 'HO2' 
# (1..., 261, ) have also been removed from the globchem

gc2.copy_species('MOH', 'MOH_A', keep_original=True)
#gc2.copy_species('HO2', ['HO2_A', 'HO2_B'], keep_original=True)


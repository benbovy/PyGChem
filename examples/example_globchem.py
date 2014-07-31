# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 18:01:16 2012

@author: bovy

Examples using the globchem module of pygchem.
See also the 'globchem' widget in the 'apps' folder.
"""

import os

import numpy as np

import pygchem.globchem as gl


ws = "/Volumes/data/03_geo/geoschem/GEOS-Chem-rundirs/4x5/geos5/standard"
fn_globchem = "globchem.dat"

# import global chemical mechanism from globchem.dat (TODO: import from KPP fmt)
gc = gl.Globchem.from_smvgear2(os.path.join(ws, fn_globchem))

# `gc` contains a dictionary of species (keys are 
# GEOS-Chem names and values are instances of the 'globchem.Species' class).
# For example, to get Species "MOH" :
gc.species['MOH']

# `gc` also contains dictionaries of kinectic and photolysis reactions (keys are
# reaction id = integer and values are instances of the 'globchem.Reaction'
# class)
# For example, get kinetic reaction with id=0 : 
gc.reac_kn[0]

# every instance of the 'globchem.Reaction' class are connected to an instance
# of 'glochem.ReactionRate', so that reaction rate parameters easily accessible
# and reaction constants are simply calculated given temperature (and other
# variables for specific reactions). For example :
reac0 = gc.reac_kn[0]
reac0.rate.ratep       # returns the reaction rate parameters
reac0.rate.ratef       # returns the function used to calculate rate constants. 
reac0.rate.get_ratek(290)    # calculate the rate constant at 290 K according to
                             # ratep and ratef
reac0.rate.get_ratek(np.arange(200,300))   # also possible to calculate rate
                                           # constants relative to value
                                           # sequences (using numpy array or
                                           # any iterable)

# find/select species or reactions by simple or advanced filtering
# step 1: define a function (or lambda) where arguments are attributes of
# Species/Reaction instances and that return a boolean
# step 2: call appropriate filter methods passing the function defined
# Examples :
spec_flt = lambda formula, status: "H" in formula and status == "A"
gc.filter_species(spec_flt)   # select species containing hydrogen atom(s) and
                              # with 'active' status

moh_flt = lambda reactants, products: "MOH" in reactants or "MOH" in products
moh_reac = gc.filter_reaction(moh_flt, getid=False)
for r in moh_reac:
    print r.format()         # print all kinetic reaction involving methanol

body3_flt = lambda flag: flag == "P"
body3_reac = gc.filter_reaction(body3_flt, getid=False)
for r in body3_reac:
    print r.format()         # print all pressure-dep 3-body reactions


# Modifications can easily be applied on the global chemical mechanism.
# Integrity of the mechanism is automatically preserved.
# Example 1: remove a species (all reactions involving the species are 
# automatically removed) :
gc.remove_species("ISOP")
isop_flt = lambda reactants, products: "ISOP" in reactants or "ISOP" in products
isop_reac = gc.filter_reaction(isop_flt)   # should return a length-0 list

# Example 2: copy or rename species (all reactions involving the species are
# updated)
gc.copy_species('MOH', ['MOHB', 'MOHO'], keep_original=False,
                name=["Methanol biofuel", "Methanol ocean em."])


# It is possible to create a global chemical mechanism from scratch...
new_gc = gl.Globchem()
my_spec = gl.Species("SPEC", name="weird", formula="SCO2", status="A",
                     atomic_mass=20.0)
new_gc.add_species(my_spec)
# ... also methods for defining and adding reactions


# TODO: export the global chemistry mechanism to the SMVGEARII (globchem.dat),
#       KPP... formats

# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2014 Benoît Bovy
# see license.txt for more details
#

"""
Functions and classes for handling GEOS-Chem datasets.

"""

import importlib


DEFAULT_BACKEND = 'iris'   # TODO: move this in the config module

_dbackend = importlib.import_module("pygchem.datafield_backends.{0}_backend"
                                    .format(DEFAULT_BACKEND))

CTMField = _dbackend.CTMField
Constraint = _dbackend.Constraint
AttributeConstraint = _dbackend.AttributeConstraint

load = _dbackend.load
load_field = _dbackend.load_field
load_fields = _dbackend.load_fields
load_raw = _dbackend.load_raw

load_callbacks = _dbackend.load_callbacks

save = _dbackend.save

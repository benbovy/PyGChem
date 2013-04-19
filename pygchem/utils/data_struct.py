# -*- coding: utf-8 -*-

# module uff
# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013 Benoît Bovy
# see license.txt for more details
# 
#
# Last modification: 04/2013

"""
Miscellaneous routine(s) and classes for data structures
"""

## {{{ http://code.activestate.com/recipes/473786/ (r1)
class AttrDict(dict):
    """
    A dictionary with attribute-style access. It maps attribute access to
    the real dictionary.
    """
    def __init__(self, init={}):
        dict.__init__(self, init)

    def __getstate__(self):
        return self.__dict__.items()

    def __setstate__(self, items):
        for key, val in items:
            self.__dict__[key] = val

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, dict.__repr__(self))

    def __setitem__(self, key, value):
        return super(AttrDict, self).__setitem__(key, value)

    def __getitem__(self, name):
        return super(AttrDict, self).__getitem__(name)

    def __delitem__(self, name):
        return super(AttrDict, self).__delitem__(name)

    __getattr__ = __getitem__
    __setattr__ = __setitem__

    def copy(self):
        ch = AttrDict(self)
        return ch

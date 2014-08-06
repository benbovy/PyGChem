# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2014 Benoit Bovy
# see license.txt for more details
#

"""
Base classes for CTM datafields.

"""


class CTMFieldBase(object):
    """
    Base class for a CTM Field.

    """

    @property
    def name(self):
        """
        The name of the data field. For CTM diagnostics, it should be
        'name_category'.
        """
        return NotImplementedError()

    @property
    def unit(self):
        """
        The unit of the data field.
        """
        return NotImplementedError()

    @property
    def attributes(self):
        """
        A dictionary of field attributes.
        """
        return NotImplementedError()

    @property
    def shape(self):
        """
        The shape of the data of this field.
        """
        return NotImplementedError()

    @property
    def ndim(self):
        """
        The number of dimensions in the data of this field.
        """
        return NotImplementedError()

    @property
    def data(self):
        """
        The :type:`numpy.ndarray` representing the (multi-dimensional) data of
        the field (will load the data into memory).
        """
        return NotImplementedError()

    @property
    def field(self):
        """
        The data field object handled by the backend library (or None if
        the backend library doesn't define any specific object for that field).
        """
        return None

    def copy(self, data=None):
        """Return a new copy of the data field."""
        return NotImplementedError()


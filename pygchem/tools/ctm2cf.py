# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013-2014 Benoit Bovy
# see license.txt for more details
#

"""
This module provides the capability to interpret CTM metadata according
to the 'NetCDF Climate and Forecast (CF) Metadata Conventions'

References:

    [CF]  NetCDF Climate and Forecast (CF) Metadata conventions, Version 1.6,
    December, 2011.

"""

import string


# unit names: CTM -> udunits2
UNITS_MAP_CTM2CF = (
    ('molec CO2', 'count'),
    ('molec', 'count'),
    ('atoms S', 'count'),
    ('atoms C', 'count'),
    ('ppbC', 'ppb'),        # prefix or suffix required (nb. of carbon atoms)
    ('kg C', 'kg'),         # prefix or suffix required (?)
    ('molC', 'mol'),        # prefix or suffix required ?     TODO: 
    ('gC', 'g'),            # prefix or suffix required ?
    ('kg S', 'kg'),
    ('kg OH', 'kg'),
    ('kg NO3', 'kg'),
    ('kg H2O2', 'kg'), 
    ('unitless', '1'),
    ('unitles', '1'),       # typo found in tracerinfo or diaginfo
    ('v/v', '1'),
    ('level', '1'),         # allowed in CF1.6 but not compatible with udunits2
    ('Fraction', '1'),
    ('fraction', '1'),
    ('ratio', '1'),
    ('factor', '1'),
    ('none', '1'),
    ('[percentage]', '%'),
    ('deg C', 'Celsius'),
    ('C', 'Celsius'),
    ('mm/da', 'mm/day'),    # typo in tracerinfo.dat 4/17/12
    ('kg/m2/', 'kg/m2'))    # ?? (tracerinfo.dat 6801 (line 1075)             


def get_cfcompliant_units(units, prefix='', suffix=''):
    """
    Get equivalent units that are compatible with the udunits2 library
    (thus CF-compliant).
    
    Parameters
    ----------
    units : string
        A string representation of the units.
    prefix : string
        Will be added at the beginning of the returned string
        (must be a valid udunits2 expression).
    suffix : string
        Will be added at the end of the returned string
        (must be a valid udunits2 expression).
    
    Returns
    -------
    A string representation of the conforming units.
    
    References
    ----------
    The udunits2 package : http://www.unidata.ucar.edu/software/udunits/
    
    Notes
    -----
    This function only relies on the table stored in :data:`UNITS_MAP_CTM2CF`.
    Therefore, the units string returned by this function is not certified to
    be compatible with udunits2.
    
    Examples
    --------
    >>> get_cfcompliant_units('molec/cm2')
    'count/cm2'
    >>> get_cfcompliant_units('v/v')
    '1'
    >>> get_cfcompliant_units('ppbC', prefix='3')
    '3ppb
    
    """
    compliant_units = units
    
    for gcunits, udunits in UNITS_MAP_CTM2CF:
        compliant_units = string.replace(compliant_units, gcunits, udunits)
    
    return prefix + compliant_units + suffix

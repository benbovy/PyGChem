# -*- coding: utf-8 -*-

# module grid
# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013 Benoit Bovy
#
# This module is partially based on the transcription of some IDL modules 
# of the GAMAP2 software (Martin Schultz, Bob Yantosca and Philippe Le Sager,
# Harvard University).
#
# see license.txt for more details
#
#
# Last modification: 06/2013

"""
This module contains routines for atmospheric profile modeling.

"""


import numpy as np


def mod_altitude(pressure, Pcoef=[-0.028389, -0.0493698, 0.485718, 0.278656,
                 -17.5703, 48.0926]):
    """
    Return altitude for given pressure.
    
    This function evaluates a polynomial at log10(pressure) values.
    
    Parameters
    ----------
    pressure : array-like
        pressure values [hPa].
    Pcoef : array-like
        coefficients of the polynomial (default values are for the US
        Standard Atmosphere).
    
    Returns
    -------
    altitude : array-like
        altitude values [km] (same shape than the pressure input array).
    
    See Also
    --------
    mod_pressure : Returns pressure for
        given altitude.
    mod_temperature : Returns air temperature for
        given altitude.
    
    Notes
    -----
    Default coefficient values represent a 5th degree polynomial which had
    been fitted to USSA data from 0-100 km. Accuracy is on the order of 1% for
    0-100 km and 0.5% below 30 km. This function, with default values, may thus
    produce bad results with pressure less than about 3e-4 hPa. 
    
    Examples
    --------
    >>> mod_altitude([1000, 800, 600])
    array([ 0.1065092 ,  1.95627858,  4.2060627 ])
    
    """
    pressure = np.asarray(pressure)
    altitude = np.polyval(Pcoef, np.log10(pressure.flatten()))
    return altitude.reshape(pressure.shape)


def mod_pressure(altitude, Zcoef=[1.94170e-9, -5.14580e-7, 4.57018e-5,
                 -1.55620e-3, -4.61994e-2, 2.99955]):
    """
    Return pressure for given altitude.
    
    This function evaluates a polynomial at altitudes values.
    
    Parameters
    ----------
    altitude : array-like
        altitude values [km].
    Zcoef : array-like
        coefficients of the polynomial (default values are for the US
        Standard Atmosphere).
    
    Returns
    -------
    pressure : array-like
        pressure values [hPa] (same shape than the altitude input array).
    
    See Also
    --------
    mod_altitude : Returns altitude for
        given pressure.
    mod_temperature : Returns air temperature for
        given altitude.
    
    Notes
    -----
    Default coefficient values represent a 5th degree polynomial which had
    been fitted to USSA data from 0-100 km. Accuracy is on the order of 1% for
    0-100 km and 0.5% below 30 km. This function, with default values, may thus
    produce bad results with altitude > 100 km.
    
    Examples
    --------
    >>> mod_pressure([0, 10, 20])
    array([ 998.96437334,  264.658697  ,   55.28114631])
    
    """
    altitude = np.asarray(altitude)
    pressure = np.power(10, np.polyval(Zcoef, altitude.flatten()))
    return pressure.reshape(altitude.shape)


def mod_temperature(altitude, Tcoef=[-4.43960e-7, 6.57752e-5, -3.62036e-3,
                    8.75339e-2, -6.75992e-1, -5.20534, 2.88283e2]):
    """
    Return temperature for given altitude.
    
    This function evaluates a polynomial at altitudes values.
    
    Parameters
    ----------
    altitude : array-like
        altitude values [km].
    Tcoef : array-like
        coefficients of the polynomial (default values are for the US
        Standard Atmosphere).
    
    Returns
    -------
    temperature : array-like
        temperature values [K] (same shape than the altitude input array).
    
    See Also
    --------
    mod_altitude : Returns altitude for
        given pressure.
    mod_pressure : Returns pressure for
        given altitude.
    
    Notes
    -----
    Default coefficient values represent a 6th order polynomial which had
    been fitted to USSA data from 0-50 km. Accuracy is on the order of 2 K
    below 8 km, and 5 K from 8-50 km. Note that this is less than the actual
    variation in atmospheric temperatures. This function, with default values,
    may thus produce bad results with altitude > 50 km.

    This function was designed to assign a temperature value to CTM grid boxes
    in order to allow conversion from mixing ratios to concentrations and
    vice versa.
    
    Examples
    --------
    >>> mod_temperature([0, 10, 20])
    array([ 288.283  ,  226.09426,  216.8602 ])
    
    """
    altitude = np.asarray(altitude)
    temperature = np.polyval(Tcoef, altitude.flatten())
    return temperature.reshape(altitude.shape)

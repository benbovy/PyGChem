# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013-2014 Benoit Bovy
#
# This module is partially based on the transcription of the IDL modules 
# `ctm_model.pro` and `ctm_grid.pro` of the GAMAP2 software (Martin Schultz,
# Bob Yantosca and Philippe Le Sager, Harvard University).
#
# see license.txt for more details
#

"""
A high-level API for handling grids used by GEOS-Chem.

It provides default grid specification for various models and allows to easily
get the grid coordinates. User-specific grids can also be defined.

"""

import numpy as np

from pygchem.tools import atm, gridspec
from pygchem.utils.numpyutil import broadcast_1d_array


class CTMGrid(object):
    """
    CTMGrid(model_name, resolution=[5, 4], halfpolar=True, center180=True,
            hybrid=True, Nlayers=None, Ntrop=None, Psurf=1013.25, Ptop=0.01,
            description='', model_family='', **kwargs):
    
    Set-up the grid of a CTM (2)3D model.

    Parameters
    ----------
    model_name : string
        Name of the model. If it is one of the supported models,
        (see :class:`CTMGrid`.supported_models), it is better to use
        :class:`CTMGrid`.from_model or :class:`CTMGrid`.copy_from_model
        to set-up the grid with appropriate parameter values. 
    resolution : (float, float)
        Horizontal grid resolution (lon, lat) or (DI, DJ) [degrees]
        (default: (5, 4))
    halfpolar : bool
        Indicates whether polar grid boxes span half (True) or same (False)
        latitude as all other boxes (default: True)
    center180 : bool
        True if lon grid is centered at 180 degrees (default: True)
    hybrid : bool
        indicates whether the model is a sigma-pressure hybrid (True) or
        pure sigma (False) level model (default: True).
    Nlayers : int or None
        Number of vertical model layers. This number must correspond to the
        number of layers in the model output files and is used in
        conjunction with Ptop to convert sigma levels into pressure
        altitudes. Set value to None if the model has no vertical
        layer (2D) (default: None).
    Ntrop : int or None
        Number of layers in the troposphere (default: None)
    Psurf : float
        Average surface pressure [hPa] (default: 1013.15)
    Ptop : float
        Pressure at model top [hPa] (default: 0.01)
    description : string
        Model grid description
    model_family : string
        Model family (e.g., 'GEOS' for 'GEOS5')

    Other Parameters
    ----------------
    Ap, Bp : 1-d array_like
        Parameters for computing ETA coordinates of the vertical grid
        levels, if hybrid (Ap [hPa] ; Bp [unitless]).
    csig, esig : 1-d array_like
        Pre-defined sigma coordinates the centers and the bottom edges of
        the vertical grid, if pure sigma.
    
    Attributes
    ----------
    Attributes are the same than the parameters above, except `model_name`
    which becomes :attr:`model`.
    
    Examples
    --------
    TODO: 

    """

    def __init__(self, model_name, resolution=(5, 4), halfpolar=True,
                 center180=True, hybrid=True, Nlayers=None, Ntrop=None,
                 Psurf=1013.25, Ptop=0.01, description='', model_family='',
                 **kwargs):

        self.model = model_name
        self.description = description
        self.model_family = model_family
        self.resolution = resolution
        self.halfpolar = bool(halfpolar)
        self.center180 = bool(center180)
        self.hybrid = bool(hybrid)
        self.Ap = None
        self.Bp = None
        self.esig = None
        self.csig = None
        try:
            self.Nlayers = int(Nlayers)
            self.Ntrop = int(Ntrop)
        except TypeError:
            self.Nlayers = Nlayers
            self.Ntrop = Ntrop
        self.Psurf = Psurf
        self.Ptop = Ptop

        self._lonlat_edges = None
        self._lonlat_centers = None
        self._eta_edges = None
        self._eta_centers = None
        self._sigma_edges = None
        self._sigma_centers = None
        self._pressure_edges = None
        self._pressure_centers = None
        self._altitude_edges = None
        self._altitude_centers = None
        
        for k, v in kwargs.items():
            self.__setattr__(k, v)
    
    @classmethod
    def from_model(cls, model_name, **kwargs):
        """
        Define a grid using the specifications of a given model.
        
        Parameters
        ----------
        model : string
            Name the model (see :func:`get_supported_models` for available
            model names).
            Supports multiple formats (e.g., 'GEOS5', 'GEOS-5' or 'GEOS_5').
        **kwargs : string
            Parameters that override the model  or default grid
          settings (See Other Parameters below).
        
        Returns
        -------
        A :class:`CTMGrid` object.
        
        Other Parameters
        ----------------
        resolution : (float, float)
            Horizontal grid resolution (lon, lat) or (DI, DJ) [degrees]
        Psurf : float
            Average surface pressure [hPa] (default: 1013.15)
        
        Notes
        -----
        Regridded vertical models may have several valid names (e.g.,
        'GEOS5_47L' and 'GEOS5_REDUCED' refer to the same model).
        
        """
        settings = gridspec.get_model_info(model_name)
        model = settings.pop('model_name')
        for k, v in kwargs.items():
            if k in ('resolution', 'Psurf'):
                settings[k] = v
        
        return cls(model, **settings)
    
    @classmethod
    def copy_from_model(cls, model_name, reference, **kwargs):
        """
        Set-up a user-defined grid using specifications of a reference
        grid model.
        
        Parameters
        ----------
        model_name : string
            name of the user-defined grid model.
        reference : string or :class:`CTMGrid` instance
            Name of the reference model (see :func:`get_supported_models`),
            or a :class:`CTMGrid` object from which grid set-up is copied.
        **kwargs
            Any set-up parameter which will override the settings of the
            reference model (see :class:`CTMGrid` parameters).
        
        Returns
        -------
        A :class:`CTMGrid` object.
        
        """
        if isinstance(reference, cls):
            settings = reference.__dict__.copy()
            settings.pop('model')
        else:
            settings = gridspec.get_model_info(reference)
            settings.pop('model_name')
        
        settings.update(kwargs)
        settings['reference'] = reference
        
        return cls(model_name, **settings)
    
    @property
    def lonlat_edges(self):
        """
        lon/lat coordinates of the edges of the grid boxes [degrees].
        
        Returns
        -------
        lon, lat : 1-d array_like 
        """
        if self._lonlat_edges is None:
            self._compute_lonlat()
        return self._lonlat_edges

    @lonlat_edges.setter
    def lonlat_edges(self, value):
        self._lonlat_edges = value
    
    @property
    def lonlat_centers(self):
        """
        lon/lat coordinates of the centers of the grid boxes [degrees].
        
        Returns
        -------
        lon, lat : 1-d array_like
        """
        if self._lonlat_centers is None:
            self._compute_lonlat() 
        return self._lonlat_centers

    @lonlat_centers.setter
    def lonlat_centers(self, value):    
        self._lonlat_centers = value
    
    def _compute_lonlat(self):
        """
        Compute lon/lat coordinates for grid edges and centers.
        
        """
        rlon = self.resolution[0]
        rlat = self.resolution[1]
        
        Nlon = int(360. / rlon)
        Nlat = int(180. / rlat) + self.halfpolar
        
        elon = (np.arange(Nlon + 1) * rlon - np.array(180.)
                - rlon / 2. * self.center180)
        elat = (np.arange(Nlat + 1) * rlat - np.array(90.)
                - rlat / 2. * self.halfpolar)
        elat[0] = -90.
        elat[-1] = 90.

        clon = (elon - rlon / 2.)[1:]
        clat = np.arange(Nlat) * rlat - np.array(90.)
            
        if self.halfpolar:                        # Fix grid boundaries
            clat[0] = (elat[0] + elat[1]) / 2.
            clat[-1] = - clat[0]
        else:
            clat += (elat[1] - elat[0]) / 2.
            
        self._lonlat_centers = (clon, clat)
        self._lonlat_edges = (elon, elat)
        
    @property
    def eta_edges(self):
        """
        ETA Coordinates of the bottom edges of the vertical hybrid grid
        [unitless].
        
        See Also
        --------
        get_layers : if not set by user, values are calculated using this
                method (with default `Psurf` and `Ptop` values).
        """
        if self._eta_edges is None:
            self._eta_edges = self.get_layers(var='eta',
                                              pos='edges')
        return self._eta_edges

    @eta_edges.setter
    def eta_edges(self, value):    
        self._eta_edges = value
    
    @property
    def eta_centers(self):
        """
        ETA Coordinates of the centers of the vertical hybrid grid [unitless].
        
        See Also
        --------
        get_layers : if not set by user, values are calculated using this
                method (with default `Psurf` and `Ptop` values).
        """
        if self._eta_centers is None:
            self._eta_centers = self.get_layers(var='eta',
                                                pos='centers')
        return self._eta_centers

    @eta_centers.setter
    def eta_centers(self, value):    
        self._eta_centers = value
    
    @property
    def sigma_edges(self):
        """
        Sigma coordinates of the vertical grid bottom edges [unitless].
        
        See Also
        --------
        get_layers : if not set by user, values are calculated using this
                method (with default `Psurf` and `Ptop` values).
        """
        if self._sigma_edges is None:
            self._sigma_edges = self.get_layers(var='sigma',
                                                pos='edges')
        return self._sigma_edges

    @sigma_edges.setter
    def sigma_edges(self, value):    
        self._sigma_edges = value
    
    @property
    def sigma_centers(self):
        """
        Sigma coordinates of the vertical grid centers [unitless].
        
        See Also
        --------
        get_layers : if not set by user, values are calculated using this
                method (with default `Psurf` and `Ptop` values).
        """
        if self._sigma_centers is None:
            self._sigma_centers = self.get_layers(var='sigma',
                                                  pos='centers')
        return self._sigma_centers

    @sigma_centers.setter
    def sigma_centers(self, value):    
        self._sigma_centers = value
    
    @property
    def pressure_edges(self):
        """
        Air pressure at the vertical grid bottom edges [hPa].
        
        See Also
        --------
        get_layers : if not set by user, values are calculated using this
                method (with default `Psurf` and `Ptop` values).
        """
        if self._pressure_edges is None:
            self._pressure_edges = self.get_layers(var='pressure',
                                                   pos='edges') 
        return self._pressure_edges

    @pressure_edges.setter
    def pressure_edges(self, value):    
        self._pressure_edges = value
    
    @property
    def pressure_centers(self):
        """
        Air pressure at the vertical grid centers [hPa].
        
        See Also
        --------
        get_layers : if not set by user, values are calculated using this
                method (with default `Psurf` and `Ptop` values).
        """
        if self._pressure_centers is None:
            self._pressure_centers = self.get_layers(var='pressure',
                                                     pos='centers') 
        return self._pressure_centers

    @pressure_centers.setter
    def pressure_centers(self, value):    
        self._pressure_centers = value
    
    @property
    def altitude_edges(self):
        """
        Altitude at the vertical grid bottom edges.
        
        See Also
        --------
        get_layers : if not set by user, values are calculated using this
                method (with default `Psurf` and `Ptop` values).
        """
        if self._altitude_edges is None:
            self._altitude_edges = self.get_layers(var='altitude',
                                                   pos='edges') 
        return self._altitude_edges

    @altitude_edges.setter
    def altitude_edges(self, value):    
        self._altitude_edges = value
    
    @property
    def altitude_centers(self):
        """
        Altitude at the vertical grid centers.
        
        See Also
        --------
        get_layers : if not set by user, values are calculated using this
                method (with default `Psurf` and `Ptop` values).
        """
        if self._altitude_centers is None:
            self._altitude_centers = self.get_layers(var='altitude',
                                                     pos='centers') 
        return self._altitude_centers

    @altitude_centers.setter
    def altitude_centers(self, value):    
        self._altitude_centers = value
    
    def get_layers(self, var, pos='edges', Psurf=None, Ptop=None,
                   **kwargs):
        """
        Compute scalars or coordinates associated to the vertical layers.
        
        Parameters
        ----------
        var : {'pressure', 'altitude', 'eta', 'sigma'}
            Scalar or coordinate to be returned.
        pos : {'edges', 'centers'}
            Get values either at grid centers or grid bottom edges
        Psurf : None or float or 1-d array-like or 2-d array-like
            Surface air pressure(s) [hPa]. If None, the value from the `Psurf`
            attribute of the :class:`CTMGrid` instance is used.
        Ptop : None or float
            Air pressure at the top of the modeled atmosphere [hPa]. If None,
            the value from the `Ptop` attribute of the :class:`CTMGrid`
            instance is used.
        
        Returns
        -------
        1-d, 2-d or 3-d array-like (the number of dimensions equals the number
        of dimensions of Psurf + 1). The first dimension of the returned
        array correspond to the vertical grid ordered bottom-up.
        
        Units are [hPa] for pressure, [km] for altitudes and [unitless] for
        eta and sigma coordinates.
        
        Returns None if `var` is not available (e.g., sigma coordinates for
        hybrid models or eta coordinates for pure sigma models) or if the
        `Nlayers` attribute of the :class:`CTMGrid` instance is None.
        
        Other Parameters
        ----------------
        Pcoef : array-like
            Coefficients of the polynomial used to get altitude values from
            pressure values (default values are for the US Standard
            Atmosphere).
        
        See Also
        --------
        eta_edges, eta_centers : Return ETA coordinates for
            hybrid grids.
        sigma_edges, sigma_centers : Return SIGMA coordinates for
            pure sigma grids.
        pressure_edges, pressure_centers : Return pressures for both
            grid types.
        altitude_edges, altitude_centers : Return altitudes for both
            grid types.
        prof_altitudes (:mod:pygchem.tools.`atm`) : Return altitude
            for given pressure.
        
        Notes
        -----
        For pure sigma grids, sigma coordinates are given by the :attr:`esig`
        (edges) and :attr:`csig` (centers) attributes of the :class:`CTMGrid`
        instance.
        
        For both pure sigma and hybrid grids, pressures at layers edges `L` are
        calculated as follows: 
        
        .. math:: P_e(L) = A_p(L) + B_p(L) * (P_{surf} - C_p)
        
        where
        
        :math:`P_{surf}`, :math:`P_{top}`
            Air pressures at the surface and the top of the modeled atmosphere
            (:attr:`Psurf` and :attr:`Ptop` attributes of the :class:`CTMGrid`
            instance).
        :math:`A_p(L)`, :math:`Bp(L)`
            Specified in the grid set-up (`Ap` and `Bp` attributes) for hybrid
            grids, or respectively equals :math:`P_{top}` and :attr:`esig`
            attribute for pure sigma grids.
        :math:`Cp(L)`
            equals :math:`P_{top}` for pure sigma grids or equals 0 for hybrid
            grids.
        
        Pressures at grid centers are averages of pressures at grid edges:
        
        .. math:: P_c(L) = (P_e(L) + P_e(L+1)) / 2
        
        For hybrid grids, ETA coordinates of grid edges and grid centers are
        given by;
        
        .. math:: ETA_{e}(L) = (P_e(L) - P_{top}) / (P_{surf} - P_{top})
        .. math:: ETA_{c}(L) = (P_c(L) - P_{top}) / (P_{surf} - P_{top})
        
        Altitude values are given using the
        :func:`pygchem.tools.atm.prof_altitudes` function.
        
        Examples
        --------
        TODO:
        
        """
        if self.Nlayers is None:
            return None
        
        if Psurf is None:
            Psurf = self.Psurf
        if Ptop is None:
            Ptop = self.Ptop
        
        Psurf = np.asarray(Psurf)
        output_ndims = Psurf.ndim + 1
        if output_ndims > 3:
            raise ValueError("`Psurf` argument must be a float or an array"
                             " with <= 2 dimensions (or None)")
        
        # compute all variables: takes not much memory, fast 
        # and better for code reading
        SIGe = None
        SIGc = None
        ETAe = None
        ETAc = None
        if self.hybrid:
            try:
                Ap = broadcast_1d_array(self.Ap, output_ndims)
                Bp = broadcast_1d_array(self.Bp, output_ndims)
            except AttributeError:
                raise AttributeError("Impossible to compute vertical levels,"
                                     " data is missing (Ap, Bp)")
            Cp = 0.
        else:
            try:
                Bp = SIGe = broadcast_1d_array(self.esig, output_ndims)
                SIGc = broadcast_1d_array(self.csig, output_ndims)   
            except AttributeError:
                raise AttributeError("Impossible to compute vertical levels,"
                                     " data is missing (esig, csig)")
            Ap = Cp = Ptop
        
        Pe = Ap + Bp * (Psurf - Cp)
        Pc = 0.5 * (Pe[0:-1] + Pe[1:])
        
        if self.hybrid:
            ETAe = (Pe - Ptop) / (Psurf - Ptop)
            ETAc = (Pc - Ptop) / (Psurf - Ptop)
        else:
            SIGe = SIGe * np.ones_like(Psurf)
            SIGc = SIGc * np.ones_like(Psurf)
        
        Ze = atm.prof_altitude(Pe, **kwargs)
        Zc = atm.prof_altitude(Pc, **kwargs)
        
        all_vars = {'eta_edges': ETAe,
                    'eta_centers': ETAc,
                    'sigma_edges': SIGe,
                    'sigma_centers': SIGc,
                    'pressure_edges': Pe,
                    'pressure_centers': Pc,
                    'altitude_edges': Ze,
                    'altitude_centers': Zc}
        
        try:
            return all_vars["_".join((var, pos))]
        except KeyError:
            raise ValueError("Invalid variable `var` and/or `pos` {'edges',"
                             " 'centers'}")

    def __repr__(self):
        return "<CTMGrid '%s' at %s>" % (self.model, hex(id(self)))

    def __str__(self):
        halfpolar = 'halfpolar, ' * self.halfpolar
        center180 = '180 degrees centered, ' * self.center180
        res = '%dx%d' % self.resolution
        horizontal = '%s%s%s lon-lat grid' % (halfpolar, center180, res)
        if self.Nlayers is None:
            hybrid = ''
            nlayers = 'no'
        else:
            hybrid = 'hybrid' * self.hybrid + 'pure sigma' * (not self.hybrid)
            nlayers = '%d layers ' % self.Nlayers
        vertical = '%s%s vertical grid' % (nlayers, hybrid)
        return "CTM Grid '%s': %s / %s (%s)" % (self.model, horizontal,
                                                vertical, self.description)


get_supported_models = gridspec.get_supported_models
get_model_info = gridspec.get_model_info

ctmgrid_from_model = CTMGrid.from_model
ctmgrid_copy = CTMGrid.copy_from_model

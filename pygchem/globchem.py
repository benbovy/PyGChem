# -*- coding: utf-8 -*-

# module globchem
# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2012 Benoit Bovy
# see license.txt for more details
#
# Last modification: 03/2013

"""
Module related to the global chemistry mechanism (species, reactions)
used in GEOS-Chem simulations.

This module provides classes and methods for reading, writing and
handling data used as input of the GEOS-Chem chemical solvers (SMVGEAR II,
FAST-J). The data is commonly stored in the files "globchem.dat",
"jv_spec.dat", "jv_spec_aod.dat", "jv_atms.dat", "ratj.d".

TODO: show here some examples.

"""

import collections
import types
import numbers
import inspect
import copy
import pprint
import warnings
import math

import numpy as np

from pygchem.utils.custom_decorators import mixedmethod, classproperty
import pygchem.io.globchem as iogc


class Species(object):
    """
    The :class:`Species` defines individual chemical species.
    """

    def __init__(self, id, name='', formula='', status='A',
                 absorption=False, atomic_mass=0.,
                 init_concentration=[0., 0., 0., 0.],
                 mb_groups={},
                 ord_n=0):
        """
        Build a :class:`Species` object.

        Parameters
        ----------
        id : string
            A given id for the species. It should correspond to the symbol
            of the species used in GEOS-Chem (also used as index for species
            assigned to a :class:`Globchem` object), thus avoid
            creating multiple :class:`Species` objects with the same id !
        name : string
            Name of the species (example: methanol).
        formula : string
            Chemical formula of the species.
        status : string
            Species status (see :class:`Species`.valid_status for a
            more detailed description).
        absorption : bool
            Tells whether species absorbs radiation (the species
            does not necessarily photolyse) (not used by GEOS-Chem which
            uses FAST-J instead)
        atomic_mass : float
            Atomic mass, in AMU
        init_contentration : list of 4 items (floats)
            Initial concentrations of species (default vol mixing ratio) in
            stratosphere, free troposphere over land, free troposphere over
            sea and urban regions, respectively.
        mb_groups : list of floats
            Species mass balances groups: used for mass balance capability of
            chemical solvers (SMVGEARII or KPP) (not used by GEOS-Chem)
        ord_n : int
            number used for species ordering
        """
        self.id = id
        self.name = name
        self.formula = formula
        self.status = status
        self.absorption = absorption
        self.atomic_mass = atomic_mass
        self.init_concentration = init_concentration
        self.mb_groups = mb_groups
        self.ord_n = ord_n

    @classproperty
    @classmethod
    def valid_status(cls):
        """
        Valid species status (a dictionary where keys are the status
        symbols and values are the status descriptions).
        """
        v_st = {'A': 'active',
                'I': 'inactive but used',
                'D': 'dead and not used',
                'U': 'used if urban chemistry is turned on',
                'S': 'used if stratospheric chemistry is turned on',
                'T': 'used if tropospheric chemistry is turned on',
                'R': 'used if urban and tropospheric chem. are turned on',
                'H': 'used if tropospheric and stratospheric chem. are '
                     'turned on'}
        return v_st

    @property
    def id(self):
        """
        Identifiant given to the species. It should correspond to the symbol
        of the species used in GEOS-Chem (also used as index for species
        assigned to a :class:`Globchem` object). (string)
        """
        return self._id
    @id.setter
    def id(self, val):
        if not isinstance(val, basestring):
            raise TypeError("Bad type for id: expected string, got %s"
                            % type(val).__name__)
        if len(val) not in xrange(1, 11):
            raise ValueError("Length of id string cannot contain more than "
                             "11 characters (but at least 1 character)")
        self._id = val.upper()

    @property
    def name(self):
        """
        Name of the species (example: methanol). (string)
        """
        return self._name
    @name.setter
    def name(self, val):
        if not isinstance(val, basestring):
            raise TypeError("Bad type for name: expected string, got %s"
                            % type(val).__name__)
        if len(val) not in xrange(0, 40):
            raise ValueError("Length of name string cannot contain more than "
                             "40 characters")
        self._name = val

    @property
    def formula(self):
        """
        Chemical formula of the species. (string)
        """
        return self._formula
    @formula.setter
    def formula(self, val):
        if not isinstance(val, basestring):
            raise TypeError("Bad type for formula: expected string, got %s"
                            % type(val).__name__)
        if len(val) not in xrange(0, 32):
            raise ValueError("Length of name string cannot contain more than "
                             "32 characters")
        self._formula = val.upper()

    @property
    def status(self):
        """
        Species status (see :class:`Species`.valid_status for a
        more detailed description). (string)
        """
        return self._status
    @status.setter
    def status(self, val):
        if not isinstance(val, basestring):
            raise TypeError("Bad type for status: expected string, got %s"
                            % type(val).__name__)
        if val not in self.valid_status.keys():
            raise ValueError("Bad value for status: expected one of the "
                             "following capital letters:\n%s"
                             % pprint.pformat(self.valid_status))
        self._status = val

    @property
    def absorption(self):
        """
        Tells whether species absorbs radiation (the species
        does not necessarily photolyse) (not used by GEOS-Chem which
        uses FAST-J instead) (bool)
        """
        return self._absorption
    @absorption.setter
    def absorption(self, val):
        if not isinstance(val, int):    # bool is a subclass of int
            raise TypeError("Bad type for absorption: expected bool, got %s"
                            % type(val).__name__)
        self._absorption = bool(val)

    @property
    def atomic_mass(self):
        """
        Atomic mass, in AMU (real number)
        """
        return self._atomic_mass
    @atomic_mass.setter
    def atomic_mass(self, val):
        if not isinstance(val, numbers.Real):
            raise TypeError("Bad type for atomic_mass: expected real number, "
                            "got %s" % type(val).__name__)
        self._atomic_mass = float(val)

    @property
    def init_concentration(self):
        """
        Initial concentrations of species (default vol mixing ratio) in
        stratosphere, free troposphere over land, free troposphere over
        sea and urban regions, respectively. (a sequence of 4 real numbers)
        """
        return self._init_concentration
    @init_concentration.setter
    def init_concentration(self, val):
        if (not isinstance(val, collections.Iterable) or len(val) != 4
            or not all(isinstance(v, numbers.Real) for v in val)):
            raise ValueError("init_concentration must be an sequence of 4 "
                             "real numbers")
        self._init_concentration = [float(v) for v in val]

    @property
    def mb_groups(self):
        """
        Species mass balances groups: used for mass balance capability of
        chemical solvers (SMVGEARII or KPP) (not used by GEOS-Chem) (dict)
        """
        return self._mb_groups
    @mb_groups.setter
    def mb_groups(self, val):
        if (not isinstance(val, dict)
            or not all(isinstance(v, numbers.Integral) for v in val.values())):
            raise ValueError("mb_groups must be a dictionary with int values")
        self._mb_groups = val

    @property
    def ord_n(self):
        """
        Number used for species sorting (integer)
        """
        return self._ord_n
    @ord_n.setter
    def ord_n(self, val):
#        if not isinstance(val, numbers.Integral):
#            raise TypeError("Bad type for ord_n: expected integer, got %s"
#                            % type(val).__name__)
        self._ord_n = int(val)

    @classmethod
    def from_dict(cls, species_dict):
        """
        Create a :class:`Species` object from as dictionary (`species_dict`)
        """
        spec = cls(species_dict['id'])
        for k in cls.get_attributes():
            if k not in species_dict.keys():
                warnings.warn("species attribute '%s' not in the dictionary, "
                              "default value %s used"
                              % (k, str(getattr(spec, k))))
        for k, v in species_dict.items():
            setattr(spec, k, v)
        return spec

    def to_dict(self):
        """
        Return a :class:`Species` object as a python dictionary
        """
        dict_attr_names = [k for k in self.__dict__.keys() if k[0] != '_']
        attr_names = set(dict_attr_names + self.get_attributes())

        return {k: getattr(self, k) for k in attr_names}

    @classmethod
    def get_attributes(cls):
        """
        Return a list of the names of the :class:`Species` built-in attributes 
        (i.e., properties). 
        """
        return [k for k, v in cls.__dict__.items() if type(v) == property]


class ReactionRate(object):
    """
    Class related to the rate constants(s) of one or several chemical
    reactions.
    """

    def __init__(self, ratef, ratep):
        """
        Build a :class:`ReactionRate` object.
        
        Parameters
        ----------
        ratef : callable
            callable that should return the reaction rate constant(s). 
            Could be either a built-in function or a
            user-supplied function. 
            Arguments of the function should be the variables
            (e.g. temperature, concentrations...) and Keyword arguments
            should be the parameters.
        ratep : dict
            dictionary with the parameters (name: val) to assign to the
            function ratef.
        
        Built-in functions
        ------------------
        :class:`ReactionRate`.arr_ratek
            used for 'standard' reactions for which the rate is defined
            by the modified Arrhenius equation
        TODO: list other built-in functions for special reactions
        """
        self.ratef = ratef
        self.ratep = ratep

    @classmethod
    def builtin_info(cls, key='flag', val='function'):
        """
        Return a dictionary with correspondence between a reaction flag and
        the built-in function (or its name) that should be used to calculate
        the rate constants.
        
        Parameters
        ----------
        key : 'flag', 'descr', 'function' or 'fname'
            keys of the returned dictionary
        val : 'flag', 'descr', 'function' or 'fname'
            values of the returned dictionary
        """

        # more consistent to create a Globchem* class attribute (or property)
        # to store flags and link to the ratek functions (TODO: )
        # also move get_ratef_fromflag
        # *or Reaction class (better) 

        flags = ['', 'P', 'E', 'X', 'Y', 'Z', 'A', 'B', 'C', 'D',
                 'K', 'V', 'G']
        functions = [cls.arr_ratek, cls.p_ratek, cls.e_ratek,
                     cls.x_ratek, cls.y_ratek, cls.z_ratek,
                     cls.a_ratek, cls.b_ratek, cls.c_ratek,
                     cls.d_ratek, cls.k_ratek, cls.v_ratek,
                     cls.g_ratek, ]
        fnames = [f.__func__.__name__ for f in functions]
        descr = flags    # TODO: short description of flags

        info = {'flag': flags, 'function': functions, 'fname': fnames,
                'descr': descr}
        if key not in info.keys() or val not in info.keys():
            raise ValueError("invalid argument(s) key and/or val")

        return dict(zip(info[key], info[val]))

    @classmethod
    def get_ratef_fromflag(cls, flag):
        """
        Return the built-in function that should be used to calculate
        the rate constants, given the reaction flag (empty flag is valid).
        """
        try:
            return cls.builtin_info()[flag]
        except KeyError:
            return None


    def get_ratek(self, *args):
        """
        Calculate the rate constant with the function (ratef) assigned to
        the :class:`ReactionRate` object.
        
        *args are the arguments required by the ratef function (temperature and
        occasionally air pressure, air density [M] or concentration of other
        species). 
        """
        return self.ratef(*args, **self.ratep)

    @classmethod
    def arr_ratek(cls, T, arr=[0., 0., 0.]):
        """
        Calculate the rate constant with the following "generic" formula
        (re-arranged from the modified Arrhenius equation):
        
        K(T) = A * (300 / T)^B * exp(C / T)
        
        Parameters
        ----------
        T : float
            temperature [K]
        arr : list of 3 floats
            parameters [A,B,C] of the rate constant expression
        """
        return arr[0] * (300 / T) ** arr[1] * np.exp(arr[2] / T)

    @classmethod
    def p_ratek(cls, T, cM, arrlow=[0., 0., 0.], arrhigh=[0., 0., 0.],
                falloff=[0., 0., 0.]):
        """
        Calculate the rate constant for pressure dependent 3-body reaction
        (falloff reaction).
        
        Falloff reaction: rate that first-order in [M] at low pressure 
        (three-body reaction), but becomes zero-order in [M] as [M] increases.
        
        The rate constant is given by:
        K(T, P_r) = K_Inf * (P_r / (1 + P_r)) * F(T, P_r)
        with
        P_r = K_0 * [M] / K_Inf
        and
        K_0   = K(0,T)   = arr_ratek(T, arrlow)   -> low  pressure limit (LPL)  
        K_Inf = K(Inf,T) = arr_ratek(T, arrhigh)  -> high pressure limit (HPL)
        and
        F(T, P_r) is the SRI Falloff function given by :
        F_c^(1 / (1 + (log_10 P_r)**2))
        with
        F_c = FCV                              or
            = exp(-T / FCT1)                   if FCT1 != 0 or
            = exp(-T / FCT1) + exp(-FCT2 / T)  if FCT1 != 0 and FCT2 != 0
        
        Parameters
        ----------
        T : float
            temperature [K]
        cM : float
            concentration of M  [TODO: units ? molec/cm3 ?]
        arrlow : list of 3 floats
            parameters [A,B,C] for the low pressure limit (LPL)
        arrhigh : list of 3 floats
            parameters [A,B,C] for the high pressure limit (HPL)
        falloff : list of 3 floats
            parameters [FCV,FCT1,FCT2] for the SRI falloff function
        
        See Also
        --------
        arr_ratek
        """
        K0 = cls.arr_ratek(T, arr=arrlow)
        Kinf = cls.arr_ratek(T, arr=arrhigh)

        Pr = K0 * cM / Kinf

        if falloff[2] != 0.:
            Fc = math.exp(-T / falloff[1]) + math.exp(-falloff[2] / T)
        elif falloff[1] != 0.:
            Fc = math.exp(-T / falloff[1])
        else:
            Fc = falloff[0]

        F = Fc ** (1. / (1. + math.log10(Pr) ** 2))

        K = Kinf * (Pr / (1. + Pr)) * F

        return K

    @classmethod
    def e_ratek(cls, T, cM, arr=[0., 0., 0.], arrlow=[0., 0., 0.],
                arrhigh=[0., 0., 0.], falloff=[0., 0., 0.]):
        """
        Calculate the rate constant of thermally dissociating species
        (reverse equilibrium reactions).
        
        The rate constant is given by:
        K(T, K_f) = K_f / (A * exp(C / T))
        where
        K_f is the rate constant of the reaction in reverse direction
        (this reaction must be pressure dependent, i.e., K_f must be
        returned by p_ratek)
        
        Parameters
        ----------
        T : float
            temperature [K]
        cM : float
            concentration of M  [TODO: units ? molec/cm3 ?]
        arr : list of 3 floats
            parameters [A,B,C] of the rate coefficient expression 
            (see arr_ratek). B should always be equal to 0.
        arrlow : list of 3 floats
            parameters [A,B,C] for the low pressure limit (LPL) (reverse
            direction)
        arrhigh : list of 3 floats
            parameters [A,B,C] for the high pressure limit (HPL) (reverse
            direction)
        falloff : list of 3 floats
            parameters [FCV,FCT1,FCT2] for the SRI falloff function (reverse
            direction)
        """
        Kf = cls.p_ratek(T, cM, f_arrlow, f_arrhigh, f_falloff)
        return Kf / cls.arr_ratek(T, arr)

    @classmethod
    def x_ratek(cls, T, cM, arrK0=[0., 0., 0.], arrK2=[0., 0., 0.],
                arrK3=[0., 0., 0.]):
        """
        Calculate the rate constant for OH + HNO3 reaction.
        
        The rate constant is given by:
        K(T, [M]) = K_0 + K_3 * [M] / (1 + K_3 * [M] / K_2)
        with
        K_0 = A0 * exp(C0 / T)
        K_2 = A2 * exp(C2 / T)
        K_3 = A3 * exp(C3 / T)
        
        Parameters
        ----------
        T : float
            temperature [K]
        cM : float
            concentration of M  [TODO: units ? molec/cm3 ?]
        arrK0 : list of 3 floats
            parameters [A0,B0,C0] of the rate constant expression K_0 (see 
            arr_ratek). B0 should always be equal to 0.
        arrK2 : list of 3 floats
            parameters [A2,B2,C2] of the rate constant expression K_2...
        arrK3 : list of 3 floats
            parameters [A3,B3,C3] of the rate constant expression K_3...
        """
        K0 = cls.arr_ratek(T, arrK0)
        K2 = cls.arr_ratek(T, arrK2)
        K3 = cls.arr_ratek(T, arrK3)
        return K0 + K3 * M / (1. + K3 * cM / K2)

    @classmethod
    def y_ratek(cls, T, Patm, arrK0=[0., 0., 0.]):
        """
        Calculate the rate constant for OH + CO reaction.
        
        The rate constant is given by:
        K(T, Patm) = K_0 * (1 + 0.6 * Patm)
        where
        K_0 = A
        
        Parameters
        ----------
        T : float
            temperature [K]
        Patm : float
            air pressure [atm]
        arrK0 : list of 3 floats
            parameters [A,B,C] of the rate constant expression K_0 (see 
            arr_ratek). B and C should always be equal to 0.
        
        TODO: new OH + CO rate constant from JPL2006 implemented in 
        calcrate.f
        """
        return cls.arr_ratek(T, arrK0) * (1. + 0.6 * Patm)

    @classmethod
    def z_ratek(cls, T, cM, cH2O, arrK1=[0., 0., 0.], arrK2=[0., 0., 0.]):
        """
        Calculate the rate constant for HO2/NO3 + HO2 reaction.
        (dependence on water vapor)
        
        The rate constant is given by:
        K(T, [M], [H2O]) =   (K_1 + K_2 * [M]) 
                           * (1 + 1.4e-21 * [H2O] * exp(2200 / T)
        where
        K_1 = A1 * exp(C1 / T)
        K_2 = A2 * exp(C2 / T)
        
        Parameters
        ----------
        T : float
            temperature [K]
        cM : float
            concentration of M  [TODO: units ? molec/cm3 ?]
        cH2O : float
            concentration of H2O  [TODO: units ? molec/cm3 ?]
        arrK1 : list of 3 floats
            parameters [A1,B1,C1] of the rate constant expression K_1 (see 
            arr_ratek). B1 should always be equal to 0.
        arrK2 : list of 3 floats
            parameters [A2,B2,C2] of the rate constant expression K_2....
        """
        K1 = cls.arr_ratek(T, arrK1)
        K2 = cls.arr_ratek(T, arrK2)
        K = (K1 + K2 * cM) * (1. + 1.4e-21 * cH2O * math.exp(2200. / T))

    @classmethod
    def a_ratek(cls, T, cM, arrK1=[0., 0., 0.], xcarbon=0.):
        """
        Calculate the rate constant for reaction with addition branch
        of RO2 + NO (forming organic nitrate)
        
        The rate constant is given by:
        K(T, [M]) = K_1 * NY
        where
        K_1 = arr_ratek(T, arrK1)
        NY (organic nitrate yields) = fyrno3(T, xcarbon, M)
        
        Parameters
        ----------
        T : float
            temperature [K]
        cM : float
            concentration of M  [TODO: units ? molec/cm3 ?]
        arrK1 : list of 3 floats
            parameters [A,B,C] of the rate constant expression K_1 (see 
            arr_ratek).
        xcarbon : float
            number of C atoms in RO2
        """
        K1 = cls.arr_ratek(T, arrK1)
        NY = cls.fyrno3(T, xcarbon, cM)
        return K1 * NY

    @classmethod
    def b_ratek(cls, T, cM, arrK1=[0., 0., 0.], xcarbon=0):
        """
        Calculate the rate constant for reaction with abstraction branch
        of RO2 + NO (forming RO + NO2)
        
        The rate constant is given by:
        K(T, M) = K_1 * (1 - NY)
        where
        K_1 = arr_ratek(T, arrK1)
        NY (organic nitrate yields) = fyrno3(T, xcarbon, M)
        
        Parameters
        ----------
        T : float
            temperature [K]
        cM : float
            concentration of M  [TODO: units ? molec/cm3 ?]
        arrK1 : list of 3 floats
            parameters [A,B,C] of the rate constant expression K_1 (see 
            arr_ratek).
        xcarbon : int
            number of C atoms in RO2
        """
        K1 = cls.arr_ratek(T, arrK1)
        NY = cls.fyrno3(T, xcarbon, cM)
        return K1 * (1. - NY)

    @classmethod
    def c_ratek(cls, T, cO2, arrK1):
        """
        Calculate the rate constant for reaction with branch of 
        GLYX + OH/NO3 forming HO2 + 2 * CO
        
        The rate constant is given by:
        K(T, [O2]) = K_1 * ([O2] + 3.5e18) / (2 * [O2] + 3.5e18)
        where
        K_1 is given by arr_ratek
        
        Parameters
        ----------
        T : float
            temperature [K]
        cO2 : float
            concentration of O2  [TODO: units ? molec/cm3 ?]
        arrK1 : list of 3 floats
            parameters [A,B,C] of the rate constant expression K_1 (see 
            arr_ratek)
        """
        K1 = cls.arr_ratek(T, arrK1)
        return K1 * (cO2 + 3.5e18) / (2. * cO2 + 3.5e18)

    @classmethod
    def d_ratek(cls, T, c02, arrK1):
        """
        Calculate the rate constant for reaction with branch of
        GLYX + OH/NO3 forming GLCO3
        
        The rate constant is given by:
        K(T, [O2]) = K_1 * [O2] / (2 * [O2] + 3.5e18)
        where
        K_1 is given by arr_ratek
        
        Parameters
        ----------
        T : float
            temperature [K]
        cO2 : float
            concentration of O2  [TODO: units ? molec/cm3 ?]
        arrK1 : list of 3 floats
            parameters [A,B,C] of the rate constant expression K_1 (see 
            arr_ratek)
        """
        K1 = cls.arr_ratek(T, arrK1)
        return K1 * cO2 / (2. * cO2 + 3.5e18)

    @classmethod
    def k_ratek(cls,):
        """
        TODO: 
        """
        pass

    @classmethod
    def v_ratek(cls, T, arrK1, arrK2):
        """
        Calculate the rate constant for reaction with temperature
        dependent branching ratio (e.g., MCO3 + MO2).
        
        The rate constant is given by:
        K(T) = K_1 / (1 + K_2)
        where
        K_1 and K_2 are given by arr_ratek
        """
        K1 = cls.arr_ratek(T, arrK1)
        K2 = cls.arr_ratek(T, arrK2)
        return K1 / (1. + K2)

    @classmethod
    def g_ratek(cls, T, cO2, arrK1, arrK2):
        """
        Calculate the rate constant for reaction DMS + OH + O2.
        
        The rate constant is given by:
        K(T) = K_1 / (1 + K_2 * [O2])
        where
        K_1 and K_2 are given by arr_ratek
        
        Parameters
        ----------
        T : float
            temperature [K]
        cO2 : float
            concentration of O2  [TODO: units ? molec/cm3 ?]
        arrK1 : list of 3 floats
            parameters [A,B,C] of the rate constant expression K_1 (see 
            arr_ratek)
        arrK2 : list of 3 floats
            parameters [A2,B2,C2] of the rate constant expression K_2....
        """
        K1 = cls.arr_ratek(T, arrK1)
        K2 = cls.arr_ratek(T, arrK2)
        return K1 / (1. + K2 * cO2)

    # TODO: built-in function for HOC2H4O ------> HO2 + 2CH2O  (see calcrate.f)

    # TODO: built-in function for HOC2H4O --O2--> HO2 + GLYC   (see calcrate.f)

    @classmethod
    def fyrno3(cls, T, xcarbn, zdnum):
        """
        Calculate organic nitrate yields YN = RKA / (RKA + RKB) from RO2+NO
        reactions as a function of the number of carbon atoms.
        
        See fyrno3.f of GeosChem source for additional notes.
        
        Parameters
        ----------
        T : float
            temperature (K)
        xcarbn : int
            number of C atoms in RO2
        zdnum : float
            air density (molec/cm3)
        """
        y300 = 0.826
        alpha = 1.94e-22
        beta = 0.97
        xm0 = 0.
        xminf = 8.1
        xf = 0.411

        xxyn = alpha * math.exp(beta * xcarbn) * zdnum * ((300. / T) ** xm0)
        yyyn = y300 * ((300. / T) ** xminf)
        aaa = math.log10(xxyn / yyyn)
        zzyn = 1. / (1. + aaa * aaa)
        rarb = (xxyn / (1. + (xxyn / yyyn))) * (xf ** zzyn)
        fyrno3 = rarb / (1. + rarb)

        return fyrno3


class Reaction(object):
    """
    Class that defines a chemical reaction.
    """

    def __init__(self, id, reactants=[], products=[], coefs=[], rate=None,
                 status='A', flag='', comment='', ord_n=0):
        """
        Build a :class:`Reaction` object.

        Parameters
        ----------
        id : int
            a given id for the reaction. It should be unique !
        reactants : list or any sequence of :class:`Species` ids
            list of reactants
        products : list or any sequence of :class:`Species` ids
            list of products
        coefs : sequence of real numbers
            product coefficients
        rate : :class:`ReactionRate` object or dict
            used to compute the rate of the chemical reaction. If a dictionary
            (with reaction rate parameters) is given, a new
            :class:`ReactionRate` object is created (considering also the flag 
            value)
        status : string
            Reaction status (see :class:`Reaction`.valid_status for a
            more detailed description).
        flag : string
            used for special reactions (i.e., reactions that have to be treated
            specifically in the model)
        comment : string
            any comment about the reaction
        ord_n : int
            number used for reaction ordering
        """
        self.id = id
        self.reactants = reactants
        self.products = products
        self.coefs = coefs
        self.flag = flag       # set flag before rate
        self.rate = rate
        self.status = status
        self.comment = comment
        self.ord_n = ord_n

    @classproperty
    @classmethod
    def valid_status(cls):
        """
        Valid reaction status (a dictionary where keys are the status
        symbols and values are the status descriptions).
        """
        v_st = {'A': 'active and included in all chemistry sets',
                'D': 'dead and not used',
                'U': 'reaction in urban chemistry set',
                'S': 'reaction in stratospheric chemistry set',
                'T': 'reaction in tropospheric chemistry set',
                'R': 'reaction in urban and tropospheric chemistry sets',
                'H': 'reaction in tropospheric and stratospheric chemistry '
                     'sets'}
        return v_st

    @property
    def id(self):
        """
        Identifiant given to the species. It should correspond to the symbol
        of the species used in GEOS-Chem (also used as index for species
        assigned to a :class:`Globchem` object). (string)
        """
        return self._id
    @id.setter
    def id(self, val):
        if not isinstance(val, numbers.Real):
            raise TypeError("bad type for id: expected number, got %s"
                            % type(val).__name__)
        self._id = val

    @property
    def reactants(self):
        """
        List of reactants
        """
        return self._reactants
    @reactants.setter
    def reactants(self, val):
        if (not isinstance(val, collections.Iterable)
            or not all(isinstance(v, basestring) for v in val)):
            raise ValueError("reactants must be a sequence of species ids")
        self._reactants = list(val)

    @property
    def products(self):
        """
        List of products
        """
        return self._products
    @products.setter
    def products(self, val):
        if (not isinstance(val, collections.Iterable)
            or not all(isinstance(v, basestring) for v in val)):
            raise ValueError("products must be a sequence of species ids")
        self._products = list(val)

    @property
    def coefs(self):
        """
        List of product coefficents
        """
        return self._coefs
    @coefs.setter
    def coefs(self, val):
        if (not isinstance(val, collections.Iterable)
            or (len(self.products) != 0 and len(val) != len(self.products))
            or not all(isinstance(v, numbers.Real) for v in val)):
            raise ValueError("coefs must be a sequence of real numbers "
                             "with the same length than products")
        self._coefs = [float(c) for c in val]

    @property
    def rate(self):
        """
        A :class:`ReactionRate` object used to compute the rate coefficents
        of the chemical reaction.
        A dictionary with the reaction rate parameters can alternatively be
        given: in that case a new :class:`ReactionRate` object is created 
        automatically (considering also the flag property).
        """
        return self._rate
    @rate.setter
    def rate(self, val):
        if isinstance(val, (ReactionRate, types.NoneType)):
            self._rate = val
        elif isinstance(val, dict):
            ratef = ReactionRate.get_ratef_fromflag(self.flag)
            rate_obj = ReactionRate(ratef, val)
            self._rate = rate_obj
        else:
            raise TypeError("Bad type for rate: expected ReactionRate object, "
                            "got %s" % type(val).__name__)

    @property
    def status(self):
        """
        Reaction status (see :class:`Reaction`.valid_status for a
        more detailed description). (string)
        """
        return self._status
    @status.setter
    def status(self, val):
        if not isinstance(val, basestring):
            raise TypeError("Bad type for status: expected string, got %s"
                            % type(val).__name__)
        if val not in self.valid_status.keys():
            raise ValueError("Bad value for status: expected one of the "
                             "following capital letters:\n%s"
                             % pprint.pformat(self.valid_status))
        self._status = val.upper()

    @property
    def flag(self):
        """
        A flag used to identify "special" reactions
        """
        return self._flag
    @flag.setter
    def flag(self, val):
        self._flag = val

    @property
    def comment(self):
        """
        Any comment about the reaction
        """
        return self._comment
    @comment.setter
    def comment(self, val):
        self._comment = str(val)

    @property
    def ord_n(self):
        """
        Number used for reaction ordering
        """
        return self._ord_n
    @ord_n.setter
    def ord_n(self, val):
        self._ord_n = int(val)

    @classmethod
    def from_dict(cls, reac_dict):
        """
        Create a :class:`Reaction` object from as dictionary (`reac_dict`)
        """
        reac = cls(reac_dict['id'])
        for k in cls.get_attributes():
            if k not in reac_dict.keys():
                warnings.warn("reaction attribute '%s' not in the dictionary, "
                              "default value %s used"
                              % (k, str(getattr(reac, k))))
        reac.flag = reac_dict['flag']   # must be set before rate
        for k, v in reac_dict.items():
            setattr(reac, k, v)
        return reac

    def to_dict(self):
        """
        Return a :class:`Reaction` object as a python dictionary
        """
        dict_attr_names = [k for k in self.__dict__.keys() if k[0] != '_']
        attr_names = set(dict_attr_names + self.get_attributes())

        return {k: getattr(self, k) for k in attr_names}

    @classmethod
    def get_attributes(cls):
        """
        Return a list of the names of the :class:`Reaction` built-in attributes 
        (i.e., properties). 
        """
        return [k for k, v in cls.__dict__.items() if type(v) == property]

    def format(self):
        """
        Format reaction to human readable format.
        
        Return
        ------
        A formatted string
        """
        fmt = " + ".join([re for re in self.reactants])
        fmt += " => "
        fmt += " + ".join(["%1.3f %s" % cre for cre in
                           zip(self.coefs, self.products)])
        return fmt
    
    def __str__(self):
        return self.format()


class KineticReaction(Reaction):
    """
    Class that defines a kinectic reaction
    """

    def __init__(self, id, **kwargs):
        """
        Build a :class:`KineticReaction` object. See :class:`Reaction` for
        information about the parameters and kwargs.
        """
        super(KineticReaction, self).__init__(id, **kwargs)

    @classmethod
    def get_attributes(cls):
        """
        Return a list of the names of the :class:`KineticReaction` built-in
        attributes (i.e., properties). Properties of the inherited
        :class:`Reaction` are also returned.
        """
        return Reaction.get_attributes()

class PhotolysisReaction(Reaction):
    """
    Class that defines a photolysis reaction
    """

    def __init__(self, id, **kwargs):
        """
        Build a :class:`PhotolysisReaction` object. See :class:`Reaction` for
        information about the parameters and kwargs.
        """
        super(PhotolysisReaction, self).__init__(id, **kwargs)

    @classmethod
    def get_attributes(cls):
        """
        Return a list of the names of the :class:`PhotolysisReaction` built-in
        attributes (i.e., properties). Properties of the inherited
        :class:`Reaction` are also returned.
        """
        return Reaction.get_attributes()


class Globchem(object):

    def __init__(self):
        """
        Class that defines a global chemistry mechanism.
        To import info from "globchem.dat", use the from_smvgear2() function.
        
        See also the "read_globchem_dat()" function of the pygchem.io.globchem
        module.
        """

        # attributes that will store global chemistry data
        self.description = ""
        self.mb_groups = dict()
        self.species = dict()
        self.reac_kn = dict()
        self.reac_ph = dict()
        self.splinks = list()

        # attributes related to the file "globchem.dat"
        self.filedir = None
        self.filename_dat = None
        self.filename_db = None
        self.filetype = "text"
        self._fline = 1

    @property
    def description(self):
        """
        Description of (or any comment about) the global chemistry mechanism
        """
        return self._description
    @description.setter
    def description(self, val):
        self._description = str(val)

    def add_mb_group(self, group_id, status='A', value=0):
        """
        Add one mass balance group and update species
        
        Parameters
        ----------
        group_id : string
            mass blance group identifiant
        status : string
            'A' (active) or 'D' (dead)
        value : int
            default value to assign to each species
        """
        if group_id in self.mb_groups.keys():
            raise KeyError("group_id already in mass balance groups")
        if status not in ('A', 'D'):
            raise ValueError("Bad value for status")
        if not isinstance(value, int):
            raise ValueError("Bad value for value: int required")

        self.mb_groups[group_id] = status

        for s in self.species.values():
            s.mb_groups[group_id] = value

    def remove_mb_group(self, group_id):
        """
        Remove a mass balance group and update species
        
        Parameters
        ----------
        group_id : string
            mass blance group identifiant
        """
        if group_id not in self.mb_groups.keys():
            raise KeyError("no mass balance group with group_id")
        del self.mb_groups[group_id]
        for s in self.species.values():
            del s.mb_groups[group_id]

    def filter_mb_groups_bystatus(self, status):
        """
        Return a id list with mass balance groups of 'status'   
        """
        return [k for k, v in self.mb_groups.items() if v == status]

    def add_species(self, spec):
        """
        Add one or several :class:`Species` object(s)

        Parameters
        ----------
        spec : :class:`Species` object or list
            the species to add
        """
        def add_one_species(spec1):
            if not isinstance(spec1, Species):
                raise TypeError("argument is not a Species object")
            elif spec1.id in self.species.keys():
                raise KeyError("species with key %s already exists in Globchem"
                           % spec1.id)
            self.species[spec1.id] = spec1

        if hasattr(spec, "__iter__") and not isinstance(spec, dict):
            for s in spec:
                add_one_species(s)
        else:
            add_one_species(spec)

    def remove_species(self, spec_id, raise_error=True):
        """
        Remove one or several :class:`Species` object(s).
        Update chemical reactions (and prompt a warning)

        Parameters
        ----------
        spec_id : string or list of strings
            the id(s) of the species to remove
        raise_error : bool
            if True, raise a KeyError if no species has an id that corresponds
            to spec_id
        """
        def remove_one_species(spec1_id):
            try:
                del self.species[spec1_id]
            except KeyError:
                if raise_error:
                    raise KeyError("species %s doesn't exist" % spec1_id)
            # update reactions
            for rt, rs in [[self.reac_kn, "kinetic"],
                           [self.reac_ph, "photolysis"]]:
                reac_list = [k for k, v in rt.items()
                             if spec1_id in v.reactants
                             or spec1_id in v.products]
                self.remove_reaction(reac_list, rs)
                if len(reac_list) > 0:
                    warnings.warn("%d %s reaction(s) involving species '%s' "
                                  "have also been removed from the globchem"
                                  % (len(reac_list), rs, spec1_id))

        if hasattr(spec_id, "__iter__") and not isinstance(spec_id, dict):
            for s in spec_id:
                remove_one_species(s)
        else:
            remove_one_species(spec_id)

    def copy_species(self, spec_id, new_id, keep_original=True, **kwargs):
        """
        "Clone" or "Split" one species into several species. Useful for
        example to use separate tracers to track a species from different
        sources (e.g., Millet et al. 2008).
        Update chemical reactions.

        Parameters
        ----------
        spec_id : string
            id of the species to copy
        new_id : string or list of strings
            new id(s) to assign
        keep_original : bool
            keep the original species
        **kwargs
            name and value of any property of the species to assign to the new
            species (see :class:`Species`).
            If new_id is a list of species ids, then the values of these
            keyword arguments must be the lists of values to assign to each
            new species
        """
        def copy_one_species(spec1_id, new1_id, **kwargs):
            if spec1_id not in self.species.keys():
                raise KeyError("Species %s doesn't exist" % spec1_id)

            new_spec = copy.deepcopy(self.species[spec1_id])
            new_spec.id = new1_id
            for k, v in kwargs.items():
                if k not in new_spec.to_dict().keys():
                    raise KeyError("Invalid Species property : %s" % k)
                setattr(new_spec, k, v)

            self.add_species(new_spec)

            # update reactions
            for rt, rs in [[self.reac_kn, "kinetic"],
                           [self.reac_ph, "photolysis"]]:
                reac_list1 = [k for k, v in rt.items()
                              if spec1_id in v.reactants]
                reac_list2 = [k for k, v in rt.items()
                              if spec1_id in v.products]
                for k in reac_list1:
                    new_id = k + 0.5
                    self.copy_reaction(k, new_id, rtype=rs, reorder=False)
                    rt[new_id].reactants = [new1_id if re == spec1_id else re
                                            for re in rt[k].reactants]
                for k in reac_list2:
                    new_id = k + 0.5
                    self.copy_reaction(k, new_id, rtype=rs, reorder=False)
                    rt[new_id].products = [new1_id if pr == spec1_id else pr
                                           for pr in rt[k].products]
            self._reorder_reactions()

        if hasattr(new_id, "__iter__") and not isinstance(new_id, dict):
            # verify and re-arrange kwargs
            for k, listv in kwargs.items():
                if len(listv) != len(new_id):
                    raise ValueError("Invalid value for kwarg '%s'" % k)
            tp = zip(*([(k, v) for v in listv] for k, listv in kwargs.items()))
            list_kwargs = map(dict, tp)

            if len(list_kwargs) > 0:
                for s, kw in zip(new_id, list_kwargs):
                    copy_one_species(spec_id, s, **kw)
            else:
                for s in new_id:
                    copy_one_species(spec_id, s)

        else:
            copy_one_species(spec_id, new_id, **kwargs)

        if not keep_original:
            self.remove_species(spec_id)

    def filter_species(self, fexpr, getid=True):
        """
        Filtering species by one or several of their attributes.
        
        Parameters
        ----------
        fexpr : function
            function that should return the result of a logical expression
            based on one or several Species attributes. The name of the
            function arguments must match the name of the :class:`Species`
            attributes ! For simple filters, a lambda function can be used.
        getid : bool
            if True, return a list of species ids.
            Otherwise, return a list of Species objects.
        
        Examples
        --------
        >>> gc = Globchem.from_smvgear2('path/to/globchem.dat')
        >>> gc.filter_species(lambda status: status == "D" or status == "I")
        >>> gc.filter_species(lambda status: status in ("I", "D"))
        >>> test_flt = lambda formula, status: formula in ("H2O", "CO2") 
        ...                                     and status == "D"
        >>> gc.filter_species(test_flt)
        """
        str_list_attr = "".join("s.%s, " % attr
                                for attr in inspect.getargspec(fexpr).args)
        str_fexpr = "fexpr(" + str_list_attr + ")"

        if getid:
            return [id for id, s in self.species.items() if eval(str_fexpr)]
        else:
            return [s for s in self.species.values() if eval(str_fexpr)]

    def add_reaction(self, reac, reorder=True):
        """
        Add one or several :class:`Reaction` object(s)

        Parameters
        ----------
        reac : :class:`Reaction` object or list of objects
            the reaction(s) to add. Instances of
            :class:`Reaction` class or :class:`KineticReaction` class are
            appended to the reac_kn attribute, and instances of
            :class:`PhotolysisReaction` class are appended to the reac_ph
            attribute.
        reorder : bool
            if True, recalculate the reactions ids/keys so that ordering
            by id remains concistent
        """
        def add_one_reaction(reac1):
            if isinstance(reac1, (PhotolysisReaction)):
                if reac1.id in self.reac_ph.keys():
                    raise KeyError("photolysis reaction with id/key %s already"
                                   " exists in Globchem" % reac1.id)
                else:
                    self.reac_ph[reac1.id] = reac1
            elif isinstance(reac1, (Reaction, KineticReaction)):
                if reac1.id in self.reac_kn.keys():
                    raise KeyError("kinetic reaction with id/key %s already "
                                   "exists in Globchem" % reac1.id)
                else:
                    self.reac_kn[reac1.id] = reac1
            else:
                raise TypeError("argument is not a Reaction object")

        if hasattr(reac, "__iter__") and not isinstance(reac, dict):
            for s in reac:
                add_one_reaction(s)
        else:
            add_one_reaction(reac)

        if reorder:
            self._reorder_reactions()

    def remove_reaction(self, reac_id, rtype='kinectic',
                        raise_error=True, reorder=True):
        """
        Remove one or several :class:`Reaction` object(s).

        Parameters
        ----------
        reac_id : string or list of strings
            the id(s) of the reaction(s) to remove
        rtype : 'kinetic' or 'photolysis'
            type of the reaction(s) to remove
        raise_error : bool
            if True, raise a KeyError if no reaction has an id that corresponds
            to reac_id
        reorder : bool
            if True, recalculate the reactions ids/keys so that ordering
            by id remains concistent
        """
        def remove_one_reaction(reac1_id, rtype='kinectic'):
            try:
                if rtype == 'kinetic':
                    del self.reac_kn[reac1_id]
                elif rtype == 'photolysis':
                    del self.reac_ph[reac1_id]
            except KeyError:
                if raise_error:
                    raise KeyError("reaction %s doesn't exist" % reac1_id)

        if hasattr(reac_id, "__iter__") and not isinstance(reac_id, dict):
            for r in reac_id:
                remove_one_reaction(r, rtype)
        else:
            remove_one_reaction(reac_id, rtype)

        if reorder:
            self._reorder_reactions()

    def copy_reaction(self, reac_id, new_id, rtype='kinetic',
                      keep_original=True, reorder=True, **kwargs):
        """
        "Duplicate" a reaction one or several times. Used by the copy_species
        method to update reactions.

        Parameters
        ----------
        reac_id : string
            id of the reaction to copy
        new_id : string or list of strings
            new id(s) to assign
        rtype : 'kinetic' or 'photolysis'
            type of the reaction(s) to duplicate
        keep_original : bool
            keep the original reaction
        reorder : bool
            if True, recalculate the reactions ids/keys so that ordering
            by id remains concistent
        **kwargs
            name and value of any property of the reaction to assign to the new
            reaction (see :class:`Reaction`).
            If new_id is a list of reaction ids, then the values of these
            keyword arguments must be the lists of values to assign to each
            new reaction
        """
        def copy_one_reaction(reac1_id, new1_id, rt, reorder, **kwargs):
            if reac1_id not in rt.keys():
                raise KeyError("Reaction %s doesn't exist" % spec1_id)

            new_reac = copy.deepcopy(rt[reac1_id])
            new_reac.id = new1_id
            for k, v in kwargs.items():
                if k not in new_reac.to_dict().keys():
                    raise KeyError("Invalid Reaction property : %s" % k)
                setattr(new_reac, k, v)

            self.add_reaction(new_reac, reorder)

        if rtype == 'kinetic':
            rt = self.reac_kn
        elif rtype == 'photolysis':
            rt = self.reac_ph

        if hasattr(new_id, "__iter__") and not isinstance(new_id, dict):
            # verify and re-arrange kwargs
            for k, listv in kwargs.items():
                if len(listv) != len(new_id):
                    raise ValueError("Invalid value for kwarg '%s'" % k)
            tp = zip(*([(k, v) for v in listv] for k, listv in kwargs.items()))
            list_kwargs = map(dict, tp)

            if len(list_kwargs) > 0:
                for r, kw in zip(new_id, list_kwargs):
                    copy_one_reaction(reac_id, r, rt, reorder, **kw)
            else:
                for r in new_id:
                    copy_one_reaction(reac_id, r, rt, reorder)

        else:
            copy_one_reaction(reac_id, new_id, rt, reorder, **kwargs)

        if not keep_original:
            self.remove_reaction(reac_id, rtype, reorder=reorder)

    def filter_reaction(self, fexpr, rtype='kinetic', getid=True):
        """
        Filtering reactions by one or several of their attributes.
        
        Parameters
        ----------
        fexpr : function
            function that should return the result of a logical expression
            based on one or several Reaction attributes. The name of the
            function arguments must match the name of the :class:`Reaction`
            attributes ! For simple filters, a lambda function can be used.
        rtype : 'kinetic' or 'photolysis'
            type of the reactions to wich filter is applied
        getid : bool
            if True, return a list of reaction ids.
            Otherwise, return a list of Reaction objects.
        
        Examples
        --------
        >>> gc = Globchem.from_smvgear2('path/to/globchem.dat')
        >>> moh_flt = lambda reactants, products: "MOH" in reactants 
        ...                                    or "MOH" in products
        >>> gc.filter_reaction(moh_flt)
        """
        str_list_attr = "".join("s.%s, " % attr
                                for attr in inspect.getargspec(fexpr).args)
        str_fexpr = "fexpr(" + str_list_attr + ")"

        if rtype == 'kinetic':
            rt = self.reac_kn
        elif rtype == 'photolysis':
            rt = self.reac_ph

        if getid:
            return [id for id, s in rt.items() if eval(str_fexpr)]
        else:
            return [s for s in rt.values() if eval(str_fexpr)]

    def _reorder_reactions(self):
        """
        Re-ordering of kinectic and photolysis reactions by ids/keys.
        Sort reaction ids/keys (integer or float) and re-assign an integer
        id/key for each reaction.
        """
        for rt in (self.reac_kn, self.reac_ph):
            new_rt = dict()
            new_key = 0
            for old_key in sorted(rt.keys()):
                new_rt[new_key] = copy.copy(rt[old_key])
                new_rt[new_key].id = new_key
                new_key += 1
            rt.clear()
            rt.update(new_rt)
            del new_rt
        warnings.warn("kinectic and/or photolysis reactions system has been "
                      "modified. ids (and keys of reac_kn and reac_ph "
                      "attributes) have been recalculated")

    @classmethod
    def from_smvgear2(cls, filename):
        """
        Import the chemistry mechanism defined in a 'globchem.dat' file
        (SMVGEAR II syntax) into a new Globchem object
        
        Parameters
        ----------
        filename : string
            name of the (or path to the) 'globchem.dat' file
        """
        (header, mb_groups, species,
         reac_kn, reac_ph) = iogc.read_globchem_dat(filename)

        gc = cls()
        gc.filename_dat = filename
        gc.description = header
        gc.mb_groups = mb_groups
        gc.add_species((Species.from_dict(s) for s in species.values()))
        gc.add_reaction((KineticReaction.from_dict(rkn)
                         for rkn in reac_kn.values()))
        gc.add_reaction((PhotolysisReaction.from_dict(rph)
                         for rph in reac_ph.values()))

        # TODO: add photolysis reactions

        return gc

    def to_smvgear2(self, filename):

        gc.filename_dat = filename
        pass

    def to_db(self, db_name):

        gc.filename_db = db_name
        pass

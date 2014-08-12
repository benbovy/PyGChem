# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2014 Benoît Bovy
# see license.txt for more details
#

"""
Interface to the Harvard-NASA Emissions Component (HEMCO).

"""

import os

from pygchem import config
from pygchem.utils.data_structures import record_cls, RecordList
from pygchem.utils import exceptions


DEFAULT_HEMCO_SETTINGS_PATH = os.path.join(config.PACKAGE_DATA_PATH,
                                           "HEMCO_settings")

#-----------------------------------------------------------------------------
# Record-type classes
#-----------------------------------------------------------------------------

_base_properties = (
    ('name', str, None, False,
     "Descriptive field name"),
    ('filename', str, None, False,
     "(netCDF) file name where the variable is stored"),
    ('var_name', str, None, False,
     "Name of the variable in the file"),
    ('ndim', int, None, False,
     "Number of dimensions of in the data"),
    ('units', str, None, False,
     "Data units"),
    ('timeslicer', str, None, False,
     "Used to get the time slices of interest from the data "
     "(see :func:`pygchem.tools.timeutil.strp_timeslicer`)"),
    ('datafield', None, None, False,
     "The datafield object")
)

_base_required = [p[0] for p in _base_properties if p[0] != 'datafield']


EmissionBase = record_cls(
    "EmissionBase",
    "Base emission field (data and) metadata",
    _base_properties +
    (('species', str, None, False,
      "Emitted chemical species"),
     ('category', int, None, False,
      "Emssion category (fields of the same category will be assembled based "
      "upon hierarchies)"),
     ('hierarchy', int, None, False,
      "Emission hierarchy (higher hierarchy emissions will overwrite lower "
      "hierarchy emissions, if same category)"),
     ('scale_factors', None, (), False,
      "Scale factors and masks (e.g., :class:`EmissionScale` or "
      ":class:`EmissionMask` objects) to be applied to "
      "the base emission field")),
    required_fields=_base_required + ['species', 'category', 'hierarchy']
)


EmissionScale = record_cls(
    "EmissionScale",
    "Emission scale factor (data and) metadata",
    _base_properties +
    (('operator', str, '*', False,
      "Mathematical operator (multiply '*', divide '/' or square '**2')"),
     ('fid', None, None, False,
      "An integer to identify the scale factor (if None, an ID will be "
      "further set automatically)")),
    required_fields=_base_required
)


EmissionMask = record_cls(
    "EmssionMask",
    "Emssion mask (data and) metadata",
    (('operator', str, '*', True,
      "Mathematical operator (should not be other than multiply '*')"),
     ('mask_window', None, None, False,
      "An approximate window for the mask field (Lon1/Lat1/Lon2/Lat2). "
      "Lon1/Lat1 is the lower left corner, Lon2/Lat2 is the upper "
      "right corner"),
     ('mirror', bool, False, False,
      "If True, invert the mask field (1-S)"),
     ('fid', None, None, False,
      "An integer to identify the scale factor (if None, an ID will be"
      "further set automatically)")),
    required_fields=_base_required
)


EmissionExt = record_cls(
    "EmissionExt",
    "An extension of the Harvard Emissions Component (HEMCO)",
    (('name', str, None, False,
      "Descriptive extension name"),
     ('enabled', bool, True, False,
      "True if the extension is enabled"),
     ('base_emission_fields', None, (), False,
      "Base emission fields (e.g., :class:`EmissionBase` objects) "
      "used by the extension"),
     ('species', None, None, False,
      "Chemical species handled by the extension"),
     ('settings', dict, {}, False,
      "Extension settings"),
     ('eid', None, None, False,
      "An integer to identify the extension (0 must correspond to HEMCO Core)."
      "If None, an ID will be further set automatically.")),
    required_fields=['name']
)


#-----------------------------------------------------------------------------
# Re-define a few properties (for more complex objects or type combinations)
#-----------------------------------------------------------------------------

# TODO: (var_name, ndim, units) correspondence between record-type and datafield

EmissionBase.scale_factors = property(
    lambda self: self._scale_factors,
    lambda self, value: setattr(
        self, '_scale_factors',
        RecordList(value, ref_classes=(EmissionScale, EmissionMask),
                   key_attr='name')
    ),
    doc=EmissionBase.scale_factors.__doc__
)

EmissionExt.base_emission_fields = property(
    lambda self: self._base_emission_fields,
    lambda self, value: setattr(
        self, '_base_emission_fields',
        RecordList(value, ref_classes=EmissionBase, key_attr='name'),
    ),
    doc=EmissionBase.scale_factors.__doc__
)


#-----------------------------------------------------------------------------
# Main emission class (+ functions)
#-----------------------------------------------------------------------------

class Emissions(object):
    """
    Settings for the Harvard-NASA Emissions Component (HEMCO).

    Parameters
    ----------
    extensions : list of :class:`EmissionExt` objects
        HEMCO extensions (must include HEMCO Core).
        Each extension can have base emission fields attached to it and each
        base emission field can have scale factors and masks attached to it.
    description : string
        A short description of emission settings.

    """
    def __init__(self, extensions=(), description=''):
        self._extensions = RecordList(extensions, ref_classes=EmissionExt)
        self.description = str(description)
        self.name = str(self.description)

    @property
    def extensions(self):
        """
        HEMCO extensions.

        Returns
        -------
        A list (:class:`pygchem.utils.data_structures.RecordList`) of extensions
        (:class:`EmissionExt` objects), including HEMCO Core.
        """
        return self._extensions

    @property
    def base_emission_fields(self):
        """
        Base emission fields.

        Returns
        -------
        A read-only list (:class:`pygchem.utils.data_structures.RecordList`)
        of all base emission fields attached to the emission setup.
        If a field is attached to several extensions, it appears only once
        in the list.
        """
        bef = []
        for ext in self.extensions:
            bef.extend(ext.base_emission_fields)
        return RecordList(set(bef), ref_classes=EmissionBase, read_only=True)

    @property
    def scale_factors(self):
        """
        Scale factors and masks.

        Returns
        -------
        A read-only list (:class:`pygchem.utils.data_structures.RecordList`)
        of all scale factors and masks attached to the emission setup.
        If a scale factor is attached to several base emission fields, it
        appears only once in the list.
        """
        sf = []
        for bef in self.base_emission_fields:
            sf.extend(bef.scale_factors)
        return RecordList(set(sf), ref_classes=[EmissionScale, EmissionMask],
                          read_only=True)

    def get_fids(self):
        """
        Return a list of all currently defined scale factor IDs.
        """
        return [sf.fid for sf in self.scale_factors]

    def get_eids(self):
        """
        Return a list of all currently defined HEMCO extension IDs.
        """
        return [ext.eid for ext in self.extensions]

    def check_ids(self):
        """
        Check the consistency of scale factors and extensions IDs.

        Raises
        ------
        ValueError
            in case of missing id duplicate IDs.

        """
        fids = self.get_fids()
        if None in fids or len(set(fids)) != len(fids):
            raise ValueError("Missing or duplicate scale factor fid")

        eids = self.get_eids()
        if None in eids or len(set(eids)) != len(eids):
            raise ValueError("Missing or duplicate extension eid")
        if 0 not in eids:
            raise ValueError("HEMCO Core (eid=0) not found in extensions")

    def resolve_ids(self):
        """
        Automatically resolve scale factors and extensions ID conficts

        Add new generated ID if not already set, update ID if needed.

        """
        fids = self.get_fids()
        min_fid = min(i for i in fids if i is not None)
        range_fids = range(min_fid, max(fids) + len(fids) + 1)
        duplicate_fids = set([fid for fid in fids if fids.count(fid) > 1])
        available_fids = sorted(set(range_fids) - set(fids))

        inc = 0
        kept_duplicate_ids = []
        for sf in self.scale_factors:
            if sf.fid not in duplicate_fids and sf.fid is not None:
                continue
            if sf.fid not in kept_duplicate_ids and sf.fid is not None:
                kept_duplicate_ids.append(sf.fid)
                continue
            sf.fid = available_fids[inc]
            inc += 1

        eids = self.get_eids()
        min_eid = min(i for i in eids if i is not None)
        range_eids = range(min_eid, max(eids) + len(eids) + 1)
        duplicate_eids = set([eid for eid in eids if eids.count(eid) > 1])
        available_eids = sorted(set(range_eids) - set(eids))

        inc = 0
        kept_duplicate_ids = []
        for ext in self.extensions:
            eid = ext.eid
            if eid not in duplicate_eids and eid is not None:
                continue
            if eid not in kept_duplicate_ids and eid is not None:
                kept_duplicate_ids.append(eid)
                continue
            ext.eid = available_eids[inc]
            inc += 1

    @classmethod
    def load(cls, filename):
        """
        Load emission settings from a file.

        Parameters
        ----------
        filename : string
            Name of (path to) the HEMCO settings file.

        Returns
        -------
        A :class:`Emissions` object.

        """
        #cls = read_config_file(filename)
        #return cls
        pass

    @classmethod
    def load_default(cls, settings):
        """
        Load one of the HEMCO settings available PyGChem.

        Parameters
        ----------
        settings : string
            Name of the settings (name of one of the files in
            :attr:`DEFAULT_HEMCO_SETTINGS_PATH`).

        Returns
        -------
        A :class:`Emissions` object.

        """
        settings_file = os.path.join(DEFAULT_HEMCO_SETTINGS_PATH, settings)
        return cls.load(settings_file)

    def save(self, filename, resolve_id=False):
        """
        Save emission settings to an HEMCO-formatted input file given by
        `filename`.
        """
        if resolve_id:
            self.resolve_ids()
        else:
            self.check_ids()

        #write_config_file(self, filename)

    def compute_emissions(self, time, grid):
        """
        Call here the Python-wrapped FORTRAN routine to calculate emissions
        at particular time(s) and with a given grid.

        Returns
        -------
        datafield object

        """
        raise exceptions.NotYetImplementedError()

    def __str__(self):
        return "HEMCO settings: {0}".format(self.description)

    def __repr__(self):
        return repr(self)


load_emissions = Emissions.load
load_emissions_default = Emissions.load_default

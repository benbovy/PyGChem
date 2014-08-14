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
from pygchem import datafields
from pygchem.io.hemco import read_hemco
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
    ('filename', None, None, False,
     "(netCDF) file name where the variable is stored"),
    ('var_name', None, None, False,
     "Name of the variable in the file"),
    ('ndim', int, None, False,
     "Number of dimensions of in the data"),
    ('units', str, None, False,
     "Data units"),
    ('timeslicer', str, None, False,
     "Used to get the time slices of interest from the data "
     "(see :func:`pygchem.tools.timeutil.strp_timeslicer`)"),
    ('datafield', None, None, False,
     "The data field (:class:`pygchem.datafields.CTMField`) object, "
     "or None if no data field has been loaded/assigned")
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
    _base_properties +
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
     ('extension_data', None, (), False,
      "Other data fields (e.g., :class:`EmissionBase` objects) "
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

EmissionExt.extension_data = property(
    lambda self: self._extension_data,
    lambda self, value: setattr(
        self, '_extension_data',
        RecordList(value, ref_classes=EmissionBase, key_attr='name'),
    ),
    doc=EmissionBase.scale_factors.__doc__
)


#-----------------------------------------------------------------------------
# Data field related functions
#-----------------------------------------------------------------------------

def load_emission_data(emission_fields, return_data=False):
    """
    Load the data field(s) corresponding to one or more given emission
    fields (base emission, scale factors and/or masks).

    Parameters
    ----------
    emission_fields : (sequence of) emission field object(s)
        load data fields for these emission fields (:class:`EmissionBase`,
        :class:`EmissionScale` or :class:`EmissionMask` objects).
    return_data : bool
        if True, it will return the loaded data fields, in addition to
        assign it to the corresponding emission fields (:prop:`datafield`).

    Notes
    -----
    The metadata and emission fields (:prop:`var_name` and :prop:`filename`)
    is used to load the data fields.

    """
    if isinstance(emission_fields, (EmissionBase, EmissionScale, EmissionMask)):
        emission_fields = [emission_fields]
    data_fields = []

    # TODO: load data fields at once for emission fields with the same filename
    for efield in emission_fields:
        constraint = datafields.Constraint(
            cube_func=lambda cube: efield.var_name == efield.var_name
        )
        dfield = datafields.load(efield.filename, constraint)
        efield.datafield = dfield
        data_fields.append(dfield)

    if return_data:
        return data_fields


def save_emission_data(emission_fields):
    """
    Save the data field(s) of to one or more given emission fields
    to their assigned file.

    Not yet implemented.

    """
    raise exceptions.NotYetImplementedError()


def gen_time_slices(emission_field):
    """
    Generate the data field time slices of interest for an emission field.

    Parameters
    ----------
    emission_field : emission field object
        :class:`EmissionBase`, :class:`EmissionScale` or :class:`EmissionMask`
        objects. :prop:`timeslicer` is used to generate the time slices.

    Returns
    -------
    generator
        sub-data field time slices.

    Not yet implemented.

    """
    raise exceptions.NotYetImplementedError()


#-----------------------------------------------------------------------------
# Main emission class (+ functions)
#-----------------------------------------------------------------------------

class EmissionSetup(object):
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
        self._extensions = RecordList(extensions, ref_classes=EmissionExt,
                                      key_attr='name')
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
        return RecordList(set(bef), ref_classes=EmissionBase, read_only=True,
                          key_attr='name')

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
                          read_only=True, key_attr='name')

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
        A :class:`EmissionSetup` object.

        """
        (settings, base_emission_fields, scale_factors, masks,
         extensions, extension_data) = read_hemco(filename)

        # add Core extension if needed
        if not any([ext['eid'] == 0 for ext in extensions]):
            extensions.append({'eid': 0, 'name': 'Core', 'enabled': True})

        emission_scales = [EmissionScale(**sf) for sf in scale_factors]
        emission_masks = [EmissionMask(**m) for m in masks]

        bef_eid_map = []
        for bef in base_emission_fields:
            bef_eid_map.append(bef.pop('eid'))
            fids = bef.pop('fids')
            bef['scale_factors'] = [sf for sf in emission_scales
                                    if sf.fid in fids]
            bef['scale_factors'] += [m for m in emission_masks if m.fid in fids]
        emission_bases = [EmissionBase(**bef) for bef in base_emission_fields]

        ed_eid_map = []
        for ed in extension_data:
            ed_eid_map.append(ed.pop('eid'))
            fids = ed.pop('fids')
            ed['scale_factors'] = [sf for sf in emission_scales
                                   if sf.fid in fids]
            ed['scale_factors'] += [m for m in emission_masks if m.fid in fids]
        emission_extdata = [EmissionBase(**ed) for ed in extension_data]

        for ext in extensions:
            ext['base_emission_fields'] = [
                bef for idx, bef in enumerate(emission_bases)
                if bef_eid_map[idx] == ext['eid']
            ]
            ext['extension_data'] = [
                ed for idx, ed in enumerate(emission_extdata)
                if ed_eid_map[idx] == ext['eid']
            ]
        emission_extensions = [EmissionExt(**ext) for ext in extensions]

        return cls(extensions=emission_extensions,
                   description="setup imported from file '{0}'"
                               .format(filename))

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
        A :class:`EmissionSetup` object.

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


load_setup = EmissionSetup.load
load_default_setup = EmissionSetup.load_default

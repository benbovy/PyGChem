# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013-2014 Christoff Keller, Benoît Bovy
# see license.txt for more details
#

"""
Read / Write Harvard-NASA Emissions Component (HEMCO) settings files.

"""

import re
import itertools
from types import StringTypes

import numpy as np


_cfg_section_wrap = '#' * 84
_cfg_comment = '#'
_cfg_extsetting_prefix = '--> '


# specific rules to convert a given string in the expected type and vice-versa
def _read_rule_srcdim(s):
    if s not in ['xy', 'xyz']:
        raise ValueError("invalid srcDim '{0}'".format(s))
    return len(s)


def _read_rule_ext_setting(s):
    if s in ['true', 'false']:
        return bool(s.title())
    try:
        return int(s)
    except Exception:
        pass
    try:
        return float(s)
    except Exception:
        pass
    return s

_cfg_read_rules = {
    'try_none': lambda s: None if s == '-' else s,
    'int': lambda s: int(s),
    'try_bool': lambda s: bool(s.title()) if s in ['true', 'false'] else s,
    'on_off': lambda s: True if s == 'on' else False,
    'ScalIDs': lambda s: [] if s == '-' else [int(i) for i in s.split('/')],
    'srcDim': _read_rule_srcdim,
    'Species': lambda s: [] if s == '-' else [sp for sp in s.split('/')],
    'mwindow': lambda s: [int(i) for i in s.split('/')],
    'Oper': lambda s: s,
    'ext_setting_name': lambda s: s.replace('--> ', ''),
    'ext_setting_val': _read_rule_ext_setting,
}

_cfg_write_rules = {
    'try_none': lambda v: '-' if v is None else str(v),
    'int': lambda v: str(v),
    'try_bool': lambda v: str(v).lower() if v is bool else str(v),
    'on_off': lambda v: 'on' if v else 'off',
    'ScalIDs': lambda v: '/'.join((str(i) for i in v)) if len(v) else '-',
    'srcDim': lambda v: 'xyz'[:v],
    'Species': lambda v: '/'.join((str(sp) for sp in v)) if len(v) else '-',
    'mwindow': lambda v: '/'.join((str(i) for i in v)),
    'ext_setting_name': lambda v: '    --> {0}'.format(v),
    'ext_setting_val': lambda v: str(v).lower() if v is bool else str(v),
}

# field separators per section
_cfg_fields_sep = {
    'SETTINGS': ':',
    'BASE EMISSIONS': ' ',
    'SCALE FACTORS': ' ',
    'MASKS': ' ',
    'EXTENSION SWITCHES': (':', ' '),   # nested splitting (':' then ' ')
    'EXTENSION SWITCHES SETTINGS': ':',
    'EXTENSION DATA': ' ',
}

# line fields specifications as a list of (field name, conversion rule),
# per section
_cfg_fields_spec = {
    'SETTINGS': (
        ('key', None),
        ('value', 'try_bool'),
    ),
    'BASE EMISSIONS': (
        ('eid', 'int'),
        ('name', None),
        ('filename', 'try_none'),
        ('var_name', 'try_none'),
        ('timeslicer', 'try_none'),
        ('ndim', 'srcDim'),
        ('units', None),
        ('species', None),
        ('fids', 'ScalIDs'),
        ('category', 'int'),
        ('hierarchy', 'int')
    ),
    'SCALE FACTORS': (
        ('fid', 'int'),
        ('name', None),
        ('filename', 'try_none'),
        ('var_name', 'try_none'),
        ('timeslicer', 'try_none'),
        ('ndim', 'srcDim'),
        ('units', None),
        ('operator', 'int')
    ),
    'MASKS': (
        ('fid', 'int'),
        ('name', None),
        ('filename', 'try_none'),
        ('var_name', 'try_none'),
        ('timeslicer', 'try_none'),
        ('ndim', 'srcDim'),
        ('units', None),
        ('operator', 'int'),
        ('mask_window', 'mwindow')
    ),
    'EXTENSION SWITCHES': (
        ('eid', 'int'),
        ('name', None),
        ('enabled', 'on_off'),
        ('species', 'Species')
    ),
    'EXTENSION SWITCHES SETTINGS': (
        ('key', 'ext_setting_name'),
        ('value', 'ext_setting_val'),
    ),

}
_cfg_fields_spec['EXTENSION DATA'] = _cfg_fields_spec['BASE EMISSIONS']


def _get_sections_lines(filename):
    """
    Returns a dictionary containing, for each config section, a list of
    (line number, raw line content) representing all configuration lines
    to be parsed (comment lines and empty lines excluded).
    """
    with open(filename) as cfg_file:

        lines = {None: []}
        section = None
        line_number = 0

        for line in cfg_file:
            line_number += 1
            line = line.strip()
            if not line or line.startswith(_cfg_comment):
                continue
            if 'BEGIN SECTION' in line:
                section = line.replace('BEGIN SECTION ', '')
                if not section in lines.keys():
                    lines[section] = []
                continue
            if 'END SECTION' in line:
                section = None
                continue

            lines[section].append((line_number, line))

    if len(lines[None]):
        raise IOError("Error while reading the HEMCO settings file '{0}': "
                      "some settings are defined outside of 'SECTION' blocks"
                      .format(filename))

    return lines


def _parse_line(line, fields_spec, fields_sep=' '):
    """
    Parse the line `line` of the configuration file, given the field
    specifications `fields_spec` and the field split character(s) `fields_sep`.
    Returns the fields as a dictionary.
    """
    # recursive line split (handle the case of EXTENSION SWITCHES)
    if isinstance(fields_sep, StringTypes):
        fields_sep = [fields_sep]
    else:
        fields_sep = list(fields_sep)

    def recursive_split(strs, seps):
        # regex (exclude any whitespace(s), tab(s)... after separator)
        sep_re = seps.pop(0) + '[\s]*'
        splitted_strs = itertools.chain.from_iterable(
            re.split(sep_re, s) for s in strs
        )
        splitted_strs = [s.strip() for s in splitted_strs]
        if len(seps):
            return recursive_split(splitted_strs, seps)
        return splitted_strs
    fields_str = recursive_split([line], fields_sep)

    if len(fields_str) != len(fields_spec):
        raise ValueError("invalid line format ('{0}')".format(fields_str))

    fields = {}
    for value, spec in zip(fields_str, fields_spec):
        key, rule_name = spec
        if rule_name is not None:
            rule = _cfg_read_rules.get(rule_name)
            value = rule(value)
        fields[key] = value

    return fields


def _parse_section_lines(lines, section, filename,
                         alt_fields_spec=None, alt_fields_sep=None):
    """parse lines for a given section and return a list of records (dicts)."""

    records = []
    fields_spec = _cfg_fields_spec[section]
    fields_sep = _cfg_fields_sep[section]

    def try_parse(ln, l, spec, sep):
        try:
            records.append(_parse_line(l, spec, sep))
        except Exception as e:
            msg = e.args[0]
            raise IOError("Error while reading the HEMCO settings file '{0}' "
                          "at line {1}: {2}"
                          .format(filename, ln, msg))

    for (line_number, line) in lines[section]:
        try:
            try_parse(line_number, line, fields_spec, fields_sep)
        except IOError:
            # allow to handle alternative specifications and separator(s)
            # (needed for extension settings)
            if alt_fields_sep is not None and alt_fields_spec is not None:
                try_parse(line_number, line, alt_fields_spec, alt_fields_sep)

    return records


def _add_datafield(efields):
    """
    Add the 'datafield' key (with appropriate value) for all given emission
    fields.
    """
    for ef in efields:
        if not len(re.sub('[0-9.\-eE/ \t]+', '', ef['filename'])):
            vals = [float(f) for f in ef['filename'].split('/')]
            if len(vals) == 1:
                ef['datafield'] = np.array(vals[0])
            else:
                ef['datafield'] = np.array(vals)
            ef['filename'] = None
        else:
            ef['datafield'] = None


def read_hemco(filename):
    """
    Read a Harvard-NASA Emissions Component (HEMCO) settings file.

    Parameters
    ----------
    filename : string
        name of (path to) the HEMCO settings file.

    Returns
    -------
    settings
        A dictionary with global HEMCO settings.
    base_emission_fields
        A list of base emission fields (metadata as dictionaries).
    scale_factors
        A list of scale factors (metadata as dictionaries).
    masks
        A list of masks metadata (metadata as dictionaries).
    extensions
        A list of HEMCO extensions (metadata as dictionaries).
    extension_data
        A list of extension data fields (metadata as dictionaries).

    """
    cfg_lines = _get_sections_lines(filename)

    base_emission_fields = _parse_section_lines(cfg_lines, 'BASE EMISSIONS',
                                                filename)
    scale_factors = _parse_section_lines(cfg_lines, 'SCALE FACTORS', filename)
    masks = _parse_section_lines(cfg_lines, 'MASKS', filename)
    extension_data = _parse_section_lines(cfg_lines, 'EXTENSION DATA', filename)

    _add_datafield(base_emission_fields + scale_factors +
                   masks + extension_data)

    # settings
    settings_list = _parse_section_lines(cfg_lines, 'SETTINGS', filename)
    settings = dict((s['key'], s['value']) for s in settings_list)

    # extensions
    extensions = _parse_section_lines(
        cfg_lines, 'EXTENSION SWITCHES', filename,
        alt_fields_spec=_cfg_fields_spec['EXTENSION SWITCHES SETTINGS'],
        alt_fields_sep=_cfg_fields_sep['EXTENSION SWITCHES SETTINGS']
    )
    ext_setting_keys = set(
        f[0] for f in _cfg_fields_spec['EXTENSION SWITCHES SETTINGS']
    )
    last_ext = None
    keep_idx = []
    for idx, ext in enumerate(extensions):
        if set(ext.keys()) == ext_setting_keys:
            last_ext['settings'][ext['key']] = ext['value']
        else:
            keep_idx.append(idx)
            last_ext = ext
            last_ext['settings'] = dict()
    extensions = [extensions[idx] for idx in keep_idx]

    return (settings, base_emission_fields, scale_factors, masks,
            extensions, extension_data)


# TODO: write_hemcp functions

# def write_config_file(emis_setup, filename, style='HEMCO'):
#     """
#     Write emission setup into file for re-use later on.
#
#     Parameters
#     ----------
#     emis_setup : :class:`Emissions` object.
#         The emissions setup object (to be written into file).
#     filename : (string)
#         File name (full path) of HEMCO configuration file
#     style : (string)
#         file type. 'HEMCO' for a HEMCO-style file, 'ESMF' for
#         an ESMF-style file.
#
#     History
#     ----------
#     20140224 ckeller: Initial version
#
#     TODO: Add ESMF writing capability.
#     """
#
#     # open file
#     outfile = open(filename, 'w')
#     outfile.write('### HEMCO INPUT FILE ###\n')
#     outfile.write(' \n')
#
#     # write file
#     write_settings(emis_setup, outfile, style)
#     write_base_emissions(emis_setup, outfile, style)
#     write_scale_factors(emis_setup, outfile, style)
#     write_masks(emis_setup, outfile, style)
#     write_extensions(emis_setup, outfile, style)
#
#     # close file
#     outfile.write('### END OF HEMCO INPUT FILE ###\n')
#     outfile.close()
#
#
# def add_header(outfile, section):
#     outfile.write('\n')
#     outfile.write(_cfg_section_wrap + '\n')
#     outfile.write('BEGIN SECTION ' + section + '\n')
#     outfile.write(_cfg_section_wrap + '\n')
#
#
# def add_footer(outfile, section):
#     outfile.write('\n')
#     outfile.write('END SECTION ' + section + '\n')
#
#
# def field2file(field, outfile, emis_setup, style, extension=None,
#                prevFile=None, prevVar=None, prevTime=None):
#     """
#     Writes a GCField object (field) to a configuration file.
#
#     Parameters
#     ----------
#     field : :class:`GCField` object.
#         field that shall be written to file.
#     outfile :
#         Already opened file to write into.
#     emis_setup : :class:`Emissions` object.
#         The emissions setup object (to be written into file).
#     style : (string)
#         file type. 'HEMCO' for a HEMCO-style file, 'ESMF' for
#         an ESMF-style file.
#     extension : :class:`EmissionExt` object.
#         The extension the field belongs to. None for scale factors.
#     prevFile, prevVar, prevTime : (optional) file name, variable and
#          time stamp of the previously written file. If the attributes
#          of the passed field match these attributes, all file values
#          (name, variable, time, unit, dimension) are set to invalid
#          value (-). This causes HEMCO to use the same file data object
#          for all fields.
#     History
#     ----------
#     20140224 ckeller: Initial version
#     20140517 ckeller: Added previous file attributes
#     """
#
#     # base field or scale field?
#     if pyhemco.emissions.BEF_ATTR_NAME in field.attributes.keys():
#         isBase = True
#         attrs = field.attributes[pyhemco.emissions.BEF_ATTR_NAME]
#     else:
#         isBase = False
#         attrs = field.attributes[pyhemco.emissions.SF_ATTR_NAME]
#
#     # extract information to be written into file
#
#     # field name
#     fldname = str(field.name)
#
#     # source file information
#     srcfile = str(field.filename)
#     srcvar = str(field.var_name)
#     srcunit = str(field.unit)
#     srctime = str(attrs['timestamp'])
#
#     # eventually write data directly into file:
#     if srcfile == '-':
#         srcfile = '/'.join([str(i) for i in list(field.data)])
#         srctime = '-'
#         srcvar = '-'
#
#     # data dimension
#     if field.ndim == 2:
#         srcdim = 'xy'
#     elif field.ndim == 3:
#         srcdim = 'xyz'
#     else:
#         raise ValueError('Illegal source dimension in ' + fldname)
#
#     # if file information are equivalent to the previous ones (passed as
#     # arguments), set all file information to invalid values (-). HEMCO
#     # will then use the same file data object for all emission fields.
#     if srcfile == prevFile and srcvar == prevVar and srctime == prevTime:
#         srcfile = '-'
#         srcvar = '-'
#         srctime = '-'
#         srcdim = '-'
#         srcunit = '-'
#
#     # BASE FIELDS
#     if isBase:
#         fid = str(extension.eid)
#         spec = str(attrs['species'])
#         cat = str(attrs['category'])
#         hier = str(attrs['hierarchy'])
#
#         scalIDs = [str(scal.attributes[pyhemco.emissions.SF_ATTR_NAME]['fid'])
#                    for scal in attrs['scale_factors']]
#         if len(scalIDs) > 0:
#             scalIDs = '/'.join(scalIDs)
#         else:
#             scalIDs = '-'
#
#         fldstr = ' '.join(
#             [fid, fldname, srcfile, srcvar, srctime, srcdim, srcunit, spec,
#              scalIDs, cat, hier])
#
#     # SCALE FACTORS / MASKS
#     else:
#         fid = str(attrs['fid'])
#         oper = str(attrs['operator'])
#         if 'mul' in oper:
#             oper = '1'
#         elif 'div' in oper:
#             oper = '-1'
#         elif 'sqr' in oper:
#             oper = '2'
#         fldstr = ' '.join(
#             [fid, fldname, srcfile, srcvar, srctime, srcdim, srcunit, oper])
#
#         # for masks:
#         if 'mask_window' in attrs.keys():
#             data = '/'.join(str(i) for i in attrs['mask_window'])
#             fldstr = ' '.join([fldstr, data])
#
#     # write to file
#     outfile.write(fldstr + '\n')
#
#
# def write_settings(emis_setup, outfile, style):
#     """
#     Write emission setup into a configuration file.
#
#     Parameters
#     ----------
#     emis_setup : :class:`Emissions` object.
#         The emissions setup object (to be written into file).
#     outfile :
#         Already opened file to write into.
#     style : (string)
#         file type. 'HEMCO' for a HEMCO-style file, 'ESMF' for
#         an ESMF-style file.
#
#     History
#     ----------
#     20140224 ckeller: Initial version
#     """
#
#     add_header(outfile, 'SETTINGS')
#     core_ext = emis_setup.extensions.get_object(name='Core')
#
#     for k, v in core_ext.settings.items():
#         outfile.write(str(k) + ': ' + str(v) + '\n')
#
#     add_footer(outfile, 'SETTINGS')
#
#
# def write_base_emissions(emis_setup, outfile, style):
#     """
#     Write base emission information into a configuration file.
#
#     Parameters
#     ----------
#     emis_setup : :class:`Emissions` object.
#         The emissions setup object (from which information is taken from).
#     outfile :
#         Already opened file to write into.
#     style : (string)
#         file type. 'HEMCO' for a HEMCO-style file, 'ESMF' for
#         an ESMF-style file.
#
#     History
#     ----------
#     20140224 ckeller: Initial version
#     """
#
#     add_header(outfile, 'BASE EMISSIONS')
#     outfile.write(
#         '# ExtNr Name sourceFile sourceVar sourceTime SrcDim SrcUnit Species ScalIDs Cat Hier\n')
#
#     core_ext = emis_setup.extensions.get_object(name='Core')
#     prevFile = ''
#     prevVar = ''
#     prevTime = ''
#     for iField in core_ext.base_emission_fields:
#         field2file(iField, outfile, emis_setup, style, core_ext, prevFile,
#                    prevVar, prevTime)
#         prevFile = str(iField.filename)
#         prevVar = str(iField.var_name)
#         prevTime = str(
#             iField.attributes[pyhemco.emissions.BEF_ATTR_NAME]['timestamp'])
#
#     add_footer(outfile, 'BASE EMISSIONS')
#
#
# def write_scale_factors(emis_setup, outfile, style):
#     """
#     Write scale factor information into a configuration file.
#
#     Parameters
#     ----------
#     emis_setup : :class:`Emissions` object.
#         The emissions setup object (from which information is taken from).
#     outfile :
#         Already opened file to write into.
#     style : (string)
#         file type. 'HEMCO' for a HEMCO-style file, 'ESMF' for
#         an ESMF-style file.
#
#     History
#     ----------
#     20140224 ckeller: Initial version
#     """
#
#     add_header(outfile, 'SCALE FACTORS')
#     outfile.write(
#         '# ScalID Name sourceFile sourceVar sourceTime SrcDim SrcUnit Oper Scalar\n')
#
#     for iField in emis_setup.scale_factors.sorted(
#             key=lambda scal: scal.attributes[pyhemco.emissions.SF_ATTR_NAME][
#                 'fid']):
#         if not iField.is_mask():
#             field2file(iField, outfile, emis_setup, style)
#
#     add_footer(outfile, 'SCALE FACTORS')
#
#
# def write_masks(emis_setup, outfile, style):
#     """
#     Write mask information into a configuration file.
#
#     Parameters
#     ----------
#     emis_setup : :class:`Emissions` object.
#         The emissions setup object (from which information is taken from).
#     outfile :
#         Already opened file to write into.
#     style : (string)
#         file type. 'HEMCO' for a HEMCO-style file, 'ESMF' for
#         an ESMF-style file.
#
#     History
#     ----------
#     20140224 ckeller: Initial version
#     """
#
#     add_header(outfile, 'MASKS')
#     outfile.write(
#         '# ScalID Name sourceFile sourceVar sourceTime SrcDim SrcUnit Oper Lon1/Lat1/Lon2/Lat2\n')
#
#     for iField in emis_setup.scale_factors.sorted(
#             key=lambda scal: scal.attributes[pyhemco.emissions.SF_ATTR_NAME][
#                 'fid']):
#         if iField.is_mask():
#             field2file(iField, outfile, emis_setup, style)
#
#     add_footer(outfile, 'MASKS')
#
#
# def write_extensions(emis_setup, outfile, style):
#     """
#     Writes extension information into a configuration file.
#
#     Parameters
#     ----------
#     emis_setup : :class:`Emissions` object.
#         The emissions setup object (from which information is taken from).
#     outfile :
#         Already opened file to write into.
#     style : (string)
#         file type. 'HEMCO' for a HEMCO-style file, 'ESMF' for
#         an ESMF-style file.
#
#     History
#     ----------
#     20140224 ckeller: Initial version
#     """
#
#     # Extension switches
#     add_header(outfile, 'EXTENSION SWITCHES')
#     outfile.write('# ExtNr ExtName     on/off Species\n')
#     for iExt in emis_setup.extensions:
#         # skip core extension: all settings already added to sections settings.
#         if iExt.eid == 0:
#             continue
#         if iExt.enabled:
#             onoff = 'on'
#         else:
#             onoff = 'off'
#         species = '/'.join(iExt.species)
#         fldstr = ' '.join([str(iExt.eid), str(iExt.name), ' :', onoff, species])
#         outfile.write(fldstr + '\n')
#         for k, v in iExt.settings.items():
#             outfile.write(
#                 '   ' + _cfg_extsetting_prefix + str(k) + ': ' + str(v) + '\n')
#     add_footer(outfile, 'EXTENSION SWITCHES')
#
#     # Extension data
#     add_header(outfile, 'EXTENSION DATA')
#     outfile.write(
#         '# ExtNr Name sourceFile sourceVar sourceTime SrcDim SrcUnit Species ScalIDs Cat Hier\n')
#     for iExt in emis_setup.extensions:
#         # skip core extension: all fields already added to section base emissions.
#         if iExt.eid == 0:
#             continue
#         for iField in iExt.base_emission_fields:
#             field2file(iField, outfile, emis_setup, style, iExt)
#     add_footer(outfile, 'EXTENSION DATA')
#
#
# def strlist_to_fields(raw_vals, fields_spec, none_val=None):
#     """
#     Get fields from a list of strings, given fields specification.
#
#     Parameters
#     ----------
#     raw_vals : [string, string, ...]
#         Values of fields that have to be convert to the correct type.
#     fields_spec : ((string, callable), (string, callable), ...)
#         Sequence of 2-length tuples (name, type) defining the name
#         and type - or any callable that return the expected type - for each
#         field.
#     none_val : string or other
#         Identifies a None value in raw_vals.
#
#     Returns
#     -------
#     tuple of 2-length tuples
#         (field name, field value).
#
#     """
#     fields = []
#     for val, spec in zip(raw_vals, field_spec):
#         fname, ftype = spec
#         if val == none_val:
#             fields.append((fname, None))
#         else:
#             fields.append((fname, ftype(val)))
#     return tuple(fields)
#
#
# def fields_to_strlist(fields, none_str=''):
#     """
#     Set a list of strings, given fields specification.
#
#     Parameters
#     ----------
#     fields : ((string, val, callable), (string, val, callable), ...)
#         (value, formatter) for each field. Formatter is a callable for
#         converting the field value to a string.
#     none_str : string
#         None value format.
#
#     Returns
#     -------
#     list of fields values as strings.
#
#     """
#     return [ffmt(fval) if fval is not None else none_str
#             for fval, ffmt in fields]
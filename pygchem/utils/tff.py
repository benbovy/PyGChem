# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013-2014 Benoît Bovy
# see license.txt for more details
#

"""
Basic functions for reading / writing formatted (Fortran) text files.

"""

from pygchem.utils.exceptions import NotYetImplementedError


def read_fmt_string(text, field_spec):
    """
    Convert a string with a given format to a dictionary.
    
    Parameters
    ----------
    text : string
        formatted string to parse
    field_spec : sequence of tuples (name, start, length, type)
        format specification for each field (field name, start position,
        char length, data type).
    
    Example
    -------
    >>> text = "Paul      28 employee"
    >>> field_spec = (('name', 0, 9, str),
    ...           ('age', 10, 12, int),
    ...           ('status', 13, str))
    >>> read_fmt_string(text, field_spec)
    {'name': 'Paul', 'age': 28, 'status': 'employee'} 
    """
    return dict((name, function(text[i:i + di].strip()))
                for name, i, di, function in field_spec)


def read_fmt_file(filename, field_spec, skip='#', as_dict=False):
    """
    Read a formatted text file.
    
    Parameters
    ----------
    filename : string
        name or path to the text file
    field_spec : sequence of tuples (name, start, length, type)
        format specification for each field (field name, start position,
        char length, data type).
    skip : string
        character(s) used to skip a line while parsing the file.
    as_dict : bool
        whether to return a list of dictionaries (records) or a dict of
        lists (fields).
    
    Returns
    -------
    if `as_dict` is False :
        a sequence (tuple) of dictionaries with field names and field
        values for each row.
    if `as_dict` is True :
        a dictionary with the field names as keys and the field values for
        all rows as values.      
    """
    with open(filename, 'r') as textfile:
        data = (read_fmt_string(line, field_spec)
                for line in textfile if not line.startswith(skip))

        if as_dict:
            fieldnames = (f[0] for f in field_spec)
            return dict(((k, list(d[k] for d in data)) for k in fieldnames))
        else:
            return list(data)


def write_fmt_file(data, filename, field_spec):
    """
    Write data to a formatted text file.

    Parameters
    ----------
    data : sequence of dicts or dict of sequences
        the structured data to write to the file.
    filename : string
        name or path to the text file.
    field_spec : sequence of tuples (name, start, length, type)
        format specification for each field (field name, start position,
        char length, data type).

    """
    # TODO: write function
    # TODO: check string length
    raise NotYetImplementedError()
# -*- coding: utf-8 -*-

# module pygchem.utils.error_handling
# pygchem: Python interface for GEOS-Chem Chemistry Transport Model
#
# Copyright (C) 2013 Benoit Bovy, Gerrit Kuhlmann
# see license.txt for more details
#
# Last modification: 04/2013

"""
Module for reading and writing text files with a fixed format
"""

def read_fmt_string(text, fields):
    """
    Convert a string with a given format to a dictionary.
    
    Parameters
    ----------
    text : string
        formatted string to parse
    fields : iterable
        sequence used to convert the text. Each item represents a field and
        must contain the field name (i.e., a key of the output dictionary),
        the start position and length of the field in text, and a conversion
        function.
    
    Example
    -------
    >>> text = "Paul      28 employee"
    >>> fields = (('name', 0, 9, str),
    ...           ('age', 10, 12, int),
    ...           ('status', 13, str))
    >>> read_fmt_string(text, fields)
    {'name': 'Paul', 'age': 28, 'status': 'employee'} 
    """
    return dict((name, function(text[i:i + di].strip()))
                 for name, i, di, function in fields)

def read_fmt_file(filename, fields, skip='#', merge_fields=False):
    """
    Read a formatted text file.
    
    Parameters
    ----------
    filename : string
        name or path to the text file
    fields : iterable
        sequence used to parse each line of the text file
        (see 'read_fmt_string')
    skip : string
        character(s) used to skip a line while parsing the file.
    
    Returns
    -------
    if merge_fields is False :
        a sequence (tuple) of dictionaries with fields names and fields
        values for each row.
    if merge_fields is True :
        a dictionary with the fields names as keys and the fields values for
        all rows as values.      
    """
    with open(filename, 'r') as textfile:
        data = (read_fmt_string(line, fields)
                for line in textfile if not line.startswith(skip))

        if merge_fields:
            fieldnames = (f[0] for f in fields)
            return dict(((k, tuple(d[k] for d in data)) for k in fieldnames))
        else:
            return tuple(data)

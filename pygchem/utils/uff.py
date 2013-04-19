# -*- coding: utf-8 -*-

# module uff
# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2012 Gerrit Kuhlmann
# see license.txt for more details
# 
#
# Last modification: 2012

"""
Module for reading and writing unformatted Fortran files.
"""


import itertools
import struct

_FIX_ERROR = (
    'Pre- and suffix of line do not match. This can happen, if the'
    ' `endian` is incorrect.'
)

class FortranFile(file):

    def __init__(self, filename, mode='rb', endian='@'):
        """
        A class for reading and writing unformatted Fortran files. Fortran
        writes data as "lines" when using the PRINT or WRITE statements.
        Each line consists of:
            - a prefix (4 byte integer gives the size of the data)
            - the real data
            - a suffix (same as prefix).

        This class can be used to read and write these "lines", in a similar
        way as reading "real lines" in a text file. A format can be given,
        while reading or writing to pack or unpack data into a binary
        format, using the 'struct' modul from the Python standard library.

        Arguments of FortranFile:

        filename, mode
            'filename' and 'mode' of the file. Default mode is 'rb' (reading
            binary).
        endian
            byte order, size and alignment of the data in
            the file. Possible characters are e.g.: '@' native (default),
            '<' little-endian, and '>' big-endian.

        See Documentation of Python's struct module for details endians and
        format strings: http://docs.python.org/library/struct.html
        """
        self.endian = endian
        file.__init__(self, filename, mode)

    def _fix(self, fmt='i'):
        """
        Read pre- or suffix of line at current position with given
        format `fmt` (default 'i').
        """
        fmt = self.endian + fmt
        fix = self.read(struct.calcsize(fmt))
        if fix:
            return struct.unpack(fmt, fix)[0]
        else:
            raise EOFError

    def readline(self, fmt=None):
        """
        Return next unformatted "line". If format is given, unpack content,
        otherwise return byte string.
        """
        size = self._fix()

        if fmt is None:
            content = self.read(size)
        else:
            fmt = self.endian + fmt
            fmt = _replace_star(fmt, size)
            content = struct.unpack(fmt, self.read(size))

        if size != self._fix():
            raise IOError(_FIX_ERROR)

        return content

    def readlines(self):
        """
        Return list strings, each a line from the file.
        """
        return [line for line in self]

    def skipline(self):
        """
        Skip the next line and returns position and size of line. Raises IOError
        if pre- and suffix of line do not match.
        """
        position = self.tell()
        prefix = self._fix()
        self.seek(prefix, 1) # skip content
        suffix = self._fix()

        if prefix != suffix:
            raise IOError(_FIX_ERROR)

        return position, prefix

    def writeline(self, fmt, *args):
        """
        Write `line` (list of objects) with given `fmt` to file. The
        `line` will be chained if object is iterable (except for
        basestrings).
        """
        fmt = self.endian + fmt
        size = struct.calcsize(fmt)

        fix = struct.pack(self.endian + 'i', size)
        line = struct.pack(fmt, *args)

        self.write(fix)
        self.write(line)
        self.write(fix)

    def writelines(self, lines, fmt):
        """
        Write `lines` with given `format`.
        """
        if isinstance(fmt, basestring):
            fmt = [fmt] * len(lines)
        for f, line in zip(fmt, lines):
            self.writeline(f, line, self.endian)

    def __iter__(self):
        return self

    def next(self, fmt=None):
        try:
            return self.readline(fmt)
        except EOFError:
            raise StopIteration


def _replace_star(fmt, size):
    """
    Replace the `*` placeholder in a format string (fmt), so that
    struct.calcsize(fmt) is equal to the given `size` using the format
    following the placeholder.

    Raises `ValueError` if number of `*`s is larger than 1. If no `*`
    in `fmt`, returns `fmt` without checking its size!

    Example
    -------
    >>> _replace_star('ii*fi', 40)
    'ii7fi'
    """
    n_stars = fmt.count('*')

    if n_stars > 1:
        raise ValueError("More than one `*` in format (%s)." % fmt)

    if n_stars:
        i = fmt.find('*')
        s = struct.calcsize(fmt.replace(fmt[i:i + 2], ''))
        n = (size - s) / struct.calcsize(fmt[i + 1])

        fmt = fmt.replace('*', str(n))

    return fmt


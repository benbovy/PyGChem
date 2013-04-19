# -*- coding: utf-8 -*-

# module pygchem.utils.check_tools
# pygchem: Python interface for GEOS-Chem Chemistry Transport Model
#
# Copyright (C) 2012 Benoit Bovy
# see license.txt for more details
#
# Last modification: 01/2013

"""
Module that defines functions for checking various things
"""

import collections


class Vref(object):
    """
    Simple class for definition of the types expected while checking a
    sequence using `check_fromref` (`ref` argument).
    """
    def __init__(self, vtype=None, length=None, vals=None):
        """
        Create a :class:`Vref` object.
        Encapsulate these objects in `ref` argument for checking an
        iterable using `check_fromref`.

        Parameters
        ----------
        vtype : None or a `type` object or a sequence of `type` objects
            allowed type(s) for the item to be checked (any type if None)
        length : None or int or sequence of int
            allowed length(s) for the item (any length if None)
        vals : list of allowed values for item (any if None)
        """
        self.vtype = vtype
        self.length = length
        self.vals = vals


def check_fromref(check, ref, raise_error=True, traceback=[]):
    """
    Recursively check elements of a sequence (dict, list, tuple or any other...
    which may itselfs encapsulate subsequences), given a reference sequence
    that contains the allowed types or classes, the allowed (sub)sequences
    lengths and the allowed values for each element.

    Parameters
    ----------
    check : sequence(s) (lists, dicts, tuples... or any iterables)
        Sequence(s) to be checked
    ref : sequence(s) (list(s) and/or dict(s))
        Reference sequence(s). Contains the allowed types or classes,
        sequences lengths, or values of each element of the sequence
        (as :class:`Vref` objects).
    raise_error : bool
        If True, raise a detailed error if check fails.
    traceback : list
        contains the base of the objects and/or indexes for error traceback 

    Returns
    -------
    bool
        True if check passes and False if check fails.

    .. note::
        The reference sequence must have the same structure of the sequence
        to be checked (and also the expected keys if sequence is and/or
        contains dictionaries). Use exclusively lists or dictionaries for
        definition of (sub)sequences in ref !
    """
    def check_type(check, types, ipath):
        """
        check element vtype subroutine
        """
        if types is None:
            return True
        if not isinstance(types, (tuple, list)):
            types = (types,)
        if not isinstance(check, types):
            if raise_error:
                fmt = "".join(["'%s' or "
                               % t.__name__ for t in types])[0:-4]
#                try:
#                    fmt = "'%s'" % iref.types.__name__
#                except AttributeError:
#                    fmt = "'%s'" % vtype(types).__name__
                raise ValueError("bad value at %s. "
                                 "Expected %s, got '%s'"
                                 % (ipath, fmt, type(check).__name__))
            return False
        else:
            return True

    def check_length(check, length, ipath):
        """
        check sequence length subroutine
        """
        if length is None:
            return True
        if isinstance(length, int):
            length = (length,)
        if len(check) not in length:
            if raise_error:
                raise ValueError("bad value at %s. Expected length = %s, "
                                 "got length = %d"
                                 % (ipath, length, len(check)))
            return False
        else:
            return True

    def check_vals(check, vals, ipath):
        """
        check element allowed values subroutine
        """
        if vals is None:
            return True
        if (not hasattr(vals, "__iter__")
              or isinstance(vals, (basestring, str))):
            vals = (vals,)
        if check not in vals:
            if raise_error:
                fmt = "".join(["'%s' or "
                               % v for v in vals])[0:-4]
                raise ValueError("bad value at %s. "
                                 "Expected one value among %s"
                                 % (ipath, fmt))
            return False
        else:
            return True

    # init check result
    check_passes = True

    # if check is dictionary: convert to sequences
    if isinstance(check, dict):
        # check a dictionary with user defined keys but values of
        # the same type (if ref is a dictionary with a unique key 'anykey')
        if 'anykey' in ref.keys():
            lcheck = check.values()
            lref = ref.values() * len(lcheck)
            lname = ["key-param '%s'" % k for k in check.keys()]
        # otherwise, check the keys 
        else:
            for k in ref.keys():
                if k not in check.keys():
                    if raise_error:
                        raise ValueError("missing key '%s' in dictionary" % k)
                    check_passes = False
                    break
            lcheck = [check[k] for k in ref.keys()]
            lref = ref.values()
            lname = ["key-param '%s'" % k for k in ref.keys()]
    # otherwise: build sequences    
    else:
        lname = ['index %i' % i for i in xrange(0, len(check))]
        lcheck = check
        lref = ref
        if len(lref) < len(lcheck):
            lref *= len(lcheck)

    # iter over elements of check and ref sequences
    for iname, icheck, iref in zip(lname, lcheck, lref):
        # format error message
        if raise_error:
            ipath = "".join([s for s in traceback]) + iname

        # element is a sequence and not a string
        # check if element must be a sequence, check length of
        # sequence, before recursive call.
        if (isinstance(icheck, collections.Iterable)
            and not isinstance(icheck, (basestring, str))):
            if not isinstance(iref, (list, dict)):
                if raise_error:
                    raise ValueError("bad value at %s. Expected %s, "
                                     "got '%s'"
                                     % (ipath,
                                        iref.vtype,
                                        type(icheck).__name__))
                check_passes = False
                break
            if not isinstance(icheck, dict):
                for i in iref:
                    if (isinstance(i, Vref)
                      and not check_length(icheck, i.length, ipath)):
                        check_passes = False
                        break

            traceback.append('%s > ' % iname)
            check_passes *= check_fromref(icheck, iref,
                                          traceback=traceback)

        # element of check is not a sequence (or string)
        #but a sequence is expected
        elif ((not isinstance(icheck, collections.Iterable)
               or isinstance(icheck, (basestring, str)))
              and isinstance(iref, (list, dict))):
            if raise_error:
                raise ValueError("bad value at %s. Expected a sequence"
                                 % icheck)
            check_passes = False
            break

        # checking element
        else:
            # check vtype and values
            check_passes *= (check_type(icheck, iref.vtype, ipath) *
                            check_vals(icheck, iref.vals, ipath))
            # if string, check length
            if isinstance(icheck, (basestring, str)):
                if not check_length(icheck, iref.length, ipath):
                    check_passes = False
                    break

    if len(traceback) > 0:
        traceback.pop(-1)  # update traceback at end of recursion step

    return check_passes

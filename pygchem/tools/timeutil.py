# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2012-2014 Gerrit Kuhlmann, Benoît Bovy
# see license.txt for more details
#

"""
Miscellaneous routine(s) for time iterations and conversions.

"""

import datetime

from dateutil.relativedelta import relativedelta


# time reference for CTM datafields
CTM_TIME_UNIT_STR = 'hours since 1985-01-01 00:00:00'
CTM_TIME_REF_DT = datetime.datetime(1985, 1, 1)


def iter_dates(start, end, step=1):
    """
    Iterate over datetime objects from `start` till `end` with `step`
    (default 1) days.

    Example:
    >>> from datetime import datetime
    >>> for date in iter_dates(datetime(2011, 1, 1), datetime(2011, 2, 1)):
    ...     print date
    """
    current = start
    while current < end:
        yield current
        current += datetime.timedelta(days=1)


def tau2time(tau, reference=CTM_TIME_REF_DT):
    """
    Convert given hours since reference (default: 01.01.1985 00:00)
    into a datetime object.
    """
    return reference + datetime.timedelta(hours=tau)


def time2tau(time, reference=CTM_TIME_REF_DT):
    """ 
    Convert a datetime object into given hours since reference
    (default: 01.01.1985 00:00).
    """
    return (time - reference).total_seconds() / 3600.0


def strp_relativedelta(delta_t):
    """
    Parse `delta_t` (string) given the format 'YYYY-MM-DD hh:mm:ss'
    and returns a :class:`dateutil.relativedelta` object.

    """
    d, t = delta_t.split(' ')
    keys = 'years', 'months', 'days', 'hours', 'minutes', 'seconds'
    vals = [int(i) for i in d.split('-')] +  [int(i) for i in t.split(':')]

    return relativedelta(**dict(zip(keys, vals)))

# -*- coding: utf-8 -*-

# module uff
# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2012-2013 Gerrit Kuhlmann, Benoît Bovy
# see license.txt for more details
# 
#
# Last modification: 04/2013

"""
Miscellaneous routine(s) for time calculations and conversions
"""

import datetime


def iter_dates(start, end, step=1):
    """
    Iterate over datetime objects from `start` till `end` with `step`
    (default 1) days.

    Example:
    >>> from datetime import datetime
    >>> for date in iter_dates(datetime(2011,1,1), datetime(2011,2,1)):
            print date
    """
    current = start
    while current < end:
        yield current
        current += datetime.timedelta(days=1)


def tau2time(tau, reference=datetime.datetime(1985, 1, 1)):
    """
    Convert given hours since reference (default: 01.01.1985 00:00)
    into a datetime object.
    """
    return reference + datetime.timedelta(hours=tau)


def time2tau(time, reference=datetime.datetime(1985, 1, 1)):
    """ 
    Convert a datetime object into given hours since reference
    (default: 01.01.1985 00:00).
    """
    return (time - reference).total_seconds() / 3600.0
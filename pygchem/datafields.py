# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2014 Benoît Bovy
# see license.txt for more details
#

"""
GEOS-Chem datasets I/O and manipulation using one of the available
backends ('iris', 'netcdf', 'bpch').

"""

import importlib
from contextlib import contextmanager


DEFAULT_BACKEND = 'iris'   # TODO: move this in the config module

# backends in failback order
_backends = ['iris', 'netcdf', 'bpch']

_current_backend = DEFAULT_BACKEND


def _load_backend(backend_name):
    global load, load_dataset, load_datasets, load_raw, load_callbacks, save
    global _current_backend

    backend = importlib.import_module("pygchem.dataset_backends.backend_{0}"
                                      .format(backend_name))

    load = backend.load
    load_dataset = backend.load_dataset
    load_datasets = backend.load_datasets
    load_raw = backend.load_raw
    load_callbacks = backend.load_callbacks
    save = backend.save

    _current_backend = backend_name


def get_available_backends():
    """
    Return a list of all available backends for dataset handling.
    """
    return _backends


@contextmanager
def set_backend(backend_name, backend_options=None, failback=True):
    """
    Set the dataset handling backend to one of the known backends.

    Parameters
    ----------
    backend_name : string
        Name of the backend
    backend_options : dict
        Not working yet
    failback : bool
        Not working yet

    See Also
    --------
    get_available_backends
        List of all available backends.
    """
    global _current_backend

    # TODO: allow to pass backend options

    backend2restore = _current_backend
    try:
        # TODO: reload backend module if backend_name == _current_backend ?
        _load_backend(backend_name)
        _current_backend = backend_name
        yield
    except Exception as e:
        # TODO: recursive call (try another backend) if failback is True
        raise e
    finally:
        _load_backend(backend2restore)


_load_backend(DEFAULT_BACKEND)

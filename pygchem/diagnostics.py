# -*- coding: utf-8 -*-

# module diagnostics
# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2012/2013 Benoit Bovy, Gerrit Kuhlmann
# see license.txt for more details
# 
#
# Last modification: 04/2013

"""
Module for reading, writing or generating GEOS-Chem diagnostics.

TODO: show here some examples.

"""

import os
import datetime

import numpy as np

from pygchem import config
from pygchem.utils import fixedformat, check_tools, uff
from pygchem.utils.data_struct import AttrDict
from pygchem.io import bpch2


class Diagnostics(object):
    """
    Define all the diagnostics (and its metadata) that one can archive
    with GEOS-Chem. An instance of this class is commonly related to a couple
    of diaginfo.dat and tracerinfo.dat files. The class provides methods to
    read/write information from/to these files.
    
    TODO: find a better way to store metadata of each diagnostic and
    category. Access metadata fields as class attributes instead than keys of
    a dictionary. Namedtuple ? a mutable, variable length version of is more 
    appropriate since we allow user-defined fields.
    Not the most elegant solution used here :
    dict subclass (AttrDict), data access as attributes (using keys as
    attr names), but no introspection (and no auto-completion possible).
    """

    diaginfo_format = (('offset', 0, 8, int),
                       ('category', 9, 40, str),
                       ('description', 49, 9999, str))

    tracerinfo_format = (('name', 0, 8, str),
                         ('full_name', 9, 30, str),
                         ('molecular_weight', 39, 10, float),
                         ('carbon_weight', 49, 3, int),
                         ('number', 52, 9, int),
                         ('scale', 61, 10, float),
                         ('unit', 72, 40, str))

    def __init__(self, diaginfo='default', tracerinfo='default'):
        """
        Initialize a :class:`Diagnostics` instance.
        
        Parameters
        ----------
        diaginfo : string
            path to the diaginfo.dat file.
            If None, no file read (no category added)
            If 'default', the instance is created using a default diaginfo.dat 
            file - located in the data folder of the module installation path)
        tracerinfo : string
            path to the tracerinfo.dat file.
            If None, no file is read (no tracer added)
            If 'default', the instance is created using the file tracerinfo.dat
            located in the data folder...
        
        See Also
        --------
        add_diagnostic, add_category (for info about diagnostic
        metadata fields)
        """
        self.categories = dict()
        self.diagnostics = dict()

        if diaginfo is not None:
            if diaginfo == 'default':
                diaginfo = os.path.join(config.PACKAGE_DATA_PATH,
                                        "diaginfo.dat")
            self.import_diaginfo(diaginfo)

        if tracerinfo is not None:
            if tracerinfo == 'default':
                tracerinfo = os.path.join(config.PACKAGE_DATA_PATH,
                                          "tracerinfo.dat")
            self.import_tracerinfo(tracerinfo)

    def import_diaginfo(self, filename, clear=True):
        """
        Read categories from the 'diaginfo.dat'-file with given `filename`.
        If `clear` is True, all existing categories are removed.
        """
        info = fixedformat.read_fmt_file(filename, self.diaginfo_format)
        if clear:
            self.categories.clear()
        self.categories.update(((d.pop('category'), AttrDict(d))
                                for d in info))

    def import_tracerinfo(self, filename, clear=True):
        """
        Read diagnostics from the 'tracerinfo.dat'-file with given `filename`.
        If `clear` is True, all existing diagnostics are removed.
        """
        info = fixedformat.read_fmt_file(filename, self.tracerinfo_format)
        if clear:
            self.diagnostics.clear()
        self.diagnostics.update(((d['number'], AttrDict(d))
                                 for d in info))

        for k, v in self.diagnostics.items():
            if v['carbon_weight'] != 1:
                self.diagnostics[k]['hydrocarbon'] = True
                self.diagnostics[k]['molecular_weight'] = (
                                            config.C_MOLECULAR_WEIGHT)
            else:
                self.diagnostics[k]['hydrocarbon'] = False

            self.diagnostics[k]['chemical'] = bool(v['molecular_weight'])

    def export_diaginfo(self):
        """
        TODO:
        """
        pass

    def export_tracerinfo(self):
        """
        TODO:
        """
        pass

    def add_category(self, category, offset, description="", **kwargs):
        """
        Add a category to the :class:`Diagnostics` instance.
        
        Parameters
        ----------
        category : string
            Category name for CTM diagnostics. The category name 
            can be up to 40 chars long, but historically the GEOS-CHEM
            and GISS models have used an 8-character category name.
        offset : int
            Constant to add to tracer numbers in order to distinguish
            for the given diagnostic category (up to 8 digits long).
        description : string
            Descriptive comment string
        
        Notes
        -----
        **kwargs are user-defined info fields assigned to the category.
        """
        self.categories.update(
                    {str(category): {'offset': int(offset),
                                     'description': str(description)}})
        self.categories[str(category)].update(kwargs)

    def del_category(self, category, raise_error=True):
        """
        Delete the `category` from the :class:`Diagnostics` instance.
        
        If `raise_error` is True, return a ValueError if `category` doesn't
        exist. Otherwise, do nothing. 
        """
        try:
            del self.categories[category]
        except KeyError:
            if raise_error:
                raise ValueError("category %s doesn't exist" % category)

    def add_category_field(self, name, value,
                           all_categories=True, category=''):
        """
        Add a metadata field of `name` and `value` for all categories or
        only `category`.
        """
        if all_categories:
            for cat in self.categories.values():
                cat.update({name: value})
        else:
            self.categories[category].update({name: value})

    def add_diagnostic(self, number, name, full_name, unit, scale, chemical,
                       molecular_weight, hydrocarbon=False, carbon_weight=0,
                       **kwargs):
        """
        Add a diagnostic to the :class:`Diagnostics` instance.
        a diagnostic can be a chemical tracer or any variable of another
        type (e.g., a MET field).
        
        Parameters
        ----------
        number : int
            Diagnostic number (up to 9 digits)
        name : string
            Diagnostic name (up to 8 chars)
        full_name : string
            Long name of the diagnostic (up to 30 chars)
        unit : string
            diagnostic unit
        scale : float
            Standard scale factor to convert values to unit given above
        chemical : bool
            True if the diagnostic is a chemical tracer
        molecular_weight : float
            Molecular weight for chemical tracers (kg/mole)
        hydrocarbon : bool
            True if the diagnostic is an hydrocarbon tracer
        carbon_weight : int
            number of moles C/moles tracer for hydrocarbon tracers
        
        Notes
        -----
        **kwargs are user-defined info fields assigned to the diagnostic. 
        """
        if hydrocarbon:
            molecular_weight = config.C_MOLECULAR_WEIGHT

        self.diagnostics.update(
                {int(number): {'name': str(name),
                               'full_name': str(full_name)},
                               'unit': str(unit),
                               'scale': float(scale),
                               'chemical': bool(chemical),
                               'molecular_weight': float(molecular_weight),
                               'hydrocarbon': bool(hydrocarbon),
                               'carbon_weight': int(carbon_weight)})
        self.diagnostics[int(number)].update(kwargs)

    def del_diagnostic(self, number, raise_error=True):
        """
        Delete the diagnostic `number` from the :class:`Diagnostics` instance.
        
        If `raise_error` is True, return a ValueError if `number` doesn't
        exist. Otherwise, do nothing. 
        """
        try:
            del self.diagnostics[number]
        except KeyError:
            if raise_error:
                raise ValueError("diagnostic %s doesn't exist" % number)

    def add_diagnostic_field(self, name, value,
                             all_diagnostics=True, number=0):
        """
        Add a metadata field of `name` and `value` for all diagnostics or
        only the diagnostic given by `number`.
        """
        if all_diagnostics:
            for diag in self.diagnostics.values():
                diag.update({name: value})
        else:
            self.diagnostics[number].update({name: value})


class DataBlock(object):
    """
    Define a data block, i.e., values of a diagnostic at a given time 
    (or relative to a given time interval).
    """

    def __init__(self, index, category, times, modelname='GEOS5_47L',
                 center180=True, halfpolar=True,
                 origin=(1, 1, 1), resolution=(5., 4.), shape=(0, 0, 0),
                 diagnostics=None, ctm_file=None, position=None, values=[],
                 **kwargs):
        """
        Create a new data block.
        
        Parameters
        ----------
        index : int
            Tracer (or diagnostic) index
        category : string
            Diagnostic category name (up to 40 chars)
        times : (datetime.datetime, datetime.datetime)
            Start and end of the diagnostic time interval
        modelname : string
            Used to identify the GEOS or GCAP MET field type relative to the
            data block (up to 20 chars)
        center180 : bool
            True if the first grid box center is located at -180 deg Lon
        halfpolar : bool
            True if the data grid has half-sized polar boxes (e.g. GEOS), False
            otherwise (e.g. GCAP)
        origin : sequence of int
            Indices of the first grid box, in Fortran notation
        resolution : (float, float)
            Lon and Lat resolution [degrees]
        shape : (int, int, int)
            shape of the data grid block
        diagnostics : object or None
            :class:`Diagnostics` instance in which the category and number of
            the diagnostic is defined. 
        ctm_file : object or None
            :class:`CTMFile` instance to which the data block is assigned
        position : int or None
            Position in the ctm_file (if format is bunch binary)
        values : numpy array or sequence of floats
            Values to assign to the data block
        
        Additional Keyword Arguments
        ----------------------------
        number : diagnostic number (should result from category offset added to
                 tracer index) 
        name : Short name of the diagnostic
        full_name : Full name of the diagnostic
        chemical : True if diagnostic is a chemical tracer
        molecular_weight : Molecular weight for chemical tracer (kg/mole)
        hydrocarbon : True if diagnostic is a hydrocarbon tracer
        carbon_weight : number of moles C/moles tracer 
                        for hydrocarbon tracers
        scale : Standard scale factor to convert values to unit
        unit : Units of the diagnostic
        
        Note
        ----
        No need to use the additional keyword arguments above if a
        :class:`Diagnostics` instance is specified.
        """
        times_ref = check_tools.Vref(datetime.datetime, 2),
        check_tools.check_fromref(times, times_ref)
        if times[1] < times[0]:
            raise ValueError("start time must be <= end time")

        self.diagnostics = diagnostics
        self.index = int(index)
        self.category = str(category)
        self.times = times
        self.modelname = str(modelname)
        self.center180 = bool(center180)
        self.halfpolar = bool(halfpolar)
        self.origin = tuple(origin)
        self.resolution = tuple(resolution)
        self._shape = tuple(shape)
        self._scale = 1.            # set both to ensure attribute exists
        self.scale = 1.             # if 'scale' is later defined as a property.

        for k, v in kwargs.items():
            setattr(self, '_' + k, v)

        # position of data in ctm_file (will be loaded on request)
        self._position = position
        self._ctm_file = ctm_file
        self._values = np.array(values)

        # verify consistency using the diagnostics instance
        if self.diagnostics is not None:
            if not isinstance(self.diagnostics, Diagnostics):
                raise ValueError("bad value for diagnostics argument")
            if self.category not in self.diagnostics.categories.keys():
                raise ValueError("bad diagnostic category")
            offset = self.diagnostics.categories[self.category]['offset']
            self._number = offset + self.index
            if self._number not in self.diagnostics.diagnostics.keys():
                raise ValueError("bad tracer index and/or diagnostic category")
            #if self._number % offset >= 100:  #100 = offset step in GEOS-Chem
            #    raise ValueError("tracer index and category don't match")

        # properties which either depend on the Diagnostics instance or
        # on attributes given by keyword arguments when calling __init__
        def _diaggetter(prop):
            def getdiag(self):
                if self.diagnostics is None:
                    return getattr(self, "_" + prop, None)
                else:
                    diag_attr = self.diagnostics.diagnostics[self._number]
                    return diag_attr.get(prop, None)
                    # see comment below for class properties
            return getdiag

        def _diagsetter(prop):
            def setdiag(self, val):
                if self.diagnostics is None:
                    setattr(self, "_" + prop, val)
                else:
                    pass   # property of Diagnostics instance is read-only here
            return setdiag

        self._inh_prop_names = [a[0] for a in Diagnostics.tracerinfo_format]
        if self.diagnostics is not None:
            self._inh_prop_names += [k for k in
                            self.diagnostics.diagnostics[self._number].keys()
                            if k not in self._inh_prop_names]
        for prop in self._inh_prop_names:
            # impossible to create properties proper to an instance. 
            # add properties at the class level. Drawback: names of attributes
            # of specific diagnostics (in the Diagnostics instance) are added
            # at the DataBlock class level (for all instances).
            setattr(self.__class__, prop, property(fget=_diaggetter(prop),
                                                   fset=_diagsetter(prop)))

    @property
    def size(self):
        #return self._values.size   # do work only if data is loaded
        return reduce(lambda x, y: x * y, self._shape)

    @property
    def shape(self):
        if len(self._values):
            return self._values.shape
        else:
            return self._shape

    @property
    def values(self):
        if (self.size == 0 and self._ctm_file is not None
            and self._position is not None):
            self._ctm_file.seek(self._position)
            vals = np.array(self._ctm_file.readline('*f'))
            vals = vals.reshape(self.shape, order='F')
            self._values = vals * self.scale
        return self._values

    @values.setter
    def values(self, data):
        self._values = np.array(data)

    def __repr__(self):
        # create time representative
        start, end = self.times
        time_fmt = '%Y-%m-%d %H:%M'
        time_repr = start.strftime(time_fmt)
        if start != end:
            time_repr = ' '.join([time_repr, '-', end.strftime(time_fmt)])

        # create name representative
        if self.name is None or self.name == '':
            name = '%s %s' % (self.number, self.category)
        else:
            name = '%s (%s) %s' % (self.name, self.number, self.category)

        return '<Datablock: %s; %s>' % (name, time_repr)


class CTMFile(object):
    """
    Represent a CTM File, and provide reading and writing methods.
    """

    supported_filetypes = ['CTM bin 02', ]

    def __init__(self, datablocks=[], diagnostics=None,
                 title='CTM Punch File by PyGChem'):
        """
        Can be used to create a :class:`CTMFile` object. To open an
        existing CTM file, use the `CTMFile.fromfile` class method
        or the `open_file` function.

        Parameters
        ----------
        datablocks : list
            List of data blocks (i.e., :class:`DataBlock` instances) to append
            to the CTM file
        diagnostics : object or None
            :class:`Diagnostics` instance to assign to the CTM file
        title : string
            CTM file title

        """

        datablocks_ref = check_tools.Vref(DataBlock),
        check_tools.check_fromref(datablocks, datablocks_ref)

        # file attributes
        self.name = None
        self.mode = None
        self.size = None
        self._file = None
        self.endian = None

        # file content
        self.datablocks = list(datablocks)
        self.diagnostics = diagnostics
        self.title = title
        self.filetype = 'Memory'

        #if self.filetype not in self.supported_filetypes:   # move to 'save'
        #    raise ValueError("Unknown filetype: %s" % self.filetype)
        if not isinstance(self.diagnostics, Diagnostics):
            raise ValueError("bad value for diagnostics argument")

    @classmethod
    def fromfile(cls, filename, mode='rb', filetype='CTM bin 02', endian='>',
                 diagnostics=None, init_diagnostics=True, skip_values=True):
        """
        Create a :class:`CTMFile` object from a given filename.

        Parameters
        ----------
        filename : string
            Name or path to the CTM file
        mode : string
            Open file mode (default 'rb')
        filetype : string
            Type of file encoding
            (default: 'CTM bin 02' for binary punch format v2)
        diagnostics : object or None
            :class:`Diagnostics` instance to assign to the CTM file
        init_diagnostics : bool
            if True, and if diagnostics is None, a new :class:`Diagnostics`
            instance is created and assigned to the CTM file, using the
            tracerinfo.dat and diaginfo.dat files in the same directory than
            filename (if they exist) or default files. 
        endian : string
            For binary punch format (default '>')
        skip_values : bool
            If True only data block headers will be read into memory and the
            values will read on request. Default: True.
        """
        # TODO: allow to load from a list of files or using a UNIX expression
        
        if init_diagnostics and diagnostics is None:
            # create a new instance of Diagnostics
            dir_path = os.path.dirname(filename)
            if dir_path == '':
                dir_path = os.getcwd()
            diaginfo = os.path.join(dir_path, "diaginfo.dat")
            tracerinfo = os.path.join(dir_path, "tracerinfo.dat")
            if not os.path.exists(diaginfo) or not os.path.exists(tracerinfo):
                diaginfo = 'default'
                tracerinfo = 'default'
            diagnostics = Diagnostics(diaginfo=diaginfo, tracerinfo=tracerinfo)

        ctm_file = cls(diagnostics=diagnostics)

        ctm_file.name = filename
        ctm_file.mode = mode
        ctm_file.datablocks = list()

        ctm_file.size = os.path.getsize(filename)
        ctm_file.endian = endian

        # read file
        if filetype == 'CTM bin 02':
            (ctm_file._file, ctm_file.filetype, ctm_file.title,
             datablocks) = bpch2.read_bpch2(filename, mode,
                                            endian, skip_values)

        # create and append datablock objects
        for datablock in datablocks:
            index = datablock.pop('index')
            category = datablock.pop('category')
            times = datablock.pop('times')
            ctm_file.append_datablock(DataBlock(index, category, times,
                                            diagnostics=ctm_file.diagnostics,
                                            **datablock))

        return ctm_file

    @property
    def categories(self):
        """
        Return a list of categories to which the data blocks contained
        in the CTM file belong.
        """
        return sorted(set(db.category for db in self.datablocks))

    def close(self):
        """
        Close file (if open).
        """
        if self._file is not None:
            self._file.close()

    def filter(self, name=None, number=None, category=None, time=None,
               fmt='%Y-%m-%d'):
        """
        Returns list of data blocks which meet 'name', 'number', 'category' and
        'time'. For 'time' a datetime.datetime object or a formatted string
        (with `fmt`) can be used.
        """
        data = self.datablocks

        if name is not None:
            data = filter(lambda db: db.name == name, data)

        if number is not None:
            data = filter(lambda db: db.number == number, data)

        if category is not None:
            data = filter(lambda db: db.category == category, data)

        if time is not None:
            if isinstance(time, basestring):
                time = datetime.datetime.strptime(time, fmt)
            data = filter(lambda db: db.times[0] == time, data)

        return sorted(data)

    def advanced_filter(self):
        """
        TODO: advanced filter based on user-defined bool function
        """

    def append_datablock(self, datablock):
        """
        append a data block (i.e., a :class:`DataBlock` instance)
        to the CTM file.
        """
        if not isinstance(datablock, DataBlock):
            raise ValueError("Invalid datablock")
        self.datablocks.append(datablock)

    def save(self, filename, filetype='CTM bin 02', endian='>',
             overwrite=False):
        """
        Save CTM file to `filename`.
        """
        if not overwrite and os.path.exists(filename):
            return True

        if filetype == 'CTM bin 02':
            bpch2.write_bpch2(self, filename, endian=endian)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if self._file is not None:
            self._file.__exit__(type, value, traceback)

    def __iter__(self):
        for obj in self.datablocks:
            yield obj

    def __len__(self):
        return len(self.datablocks)

    def __repr__(self):
        if self._file is None:
            return "<CTM file in memory>"
        else:
            status = 'closed' if self._file.closed else 'open'
            string = "<%s CTM file '%s', mode '%s' at %s>"
            return string % (status, self.name, self.mode, hex(id(self._file)))

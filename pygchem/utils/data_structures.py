# -*- coding: utf-8 -*-

# parts of pygchem (Python interface for GEOS-Chem Chemistry Transport Model)
#
# Copyright (C) 2013-2014 Benoît Bovy
# see license.txt for more details
#

"""
Custom, generic data structures.

"""

import sys
import os
from textwrap import dedent
from keyword import iskeyword
from collections import OrderedDict
from collections import Counter
try:
    from UserList import UserList
except ImportError:
    from collections import UserList

from pygchem.utils import exceptions


class RecordList(UserList):
    """
    A list of objects of one (record-type) class with basic database-like
    lookup capabilities.

    Parameters
    ----------
    items : iterable
        Objects (records) to append to the list. All items in the iterable must
        be instances of a same class.
    key_attr : string
        Name of the class attribute that will be used to define a key for each
        item in the list (if None, auto incremented integers will be used).
    triggers_add : (callable or None, callable or None)
        Functions called before / after an item is added to the list.
        Each function must accepts one argument that must be any
        object of class `ref_class`.
    triggers_del : (callable or None, callable or None)
        Functions called before / after an item is removed from the list.
    read_only : bool
        If True, the list is read-only and will behave like an immutable
        sequence.
    ref_class : class or None
        Specify explicitly the class of items in the list
        (if None, the class will be determined from the first item in `items`).

    Notes
    -----
    `key_attr` is used for uniqueness but not for indexing.

    """

    def __init__(self, items, key_attr=None,
                 triggers_add=(None, None), triggers_del=(None, None),
                 read_only=False, ref_class=None):

        super(RecordList, self).__init__([])

        self.key_attr = key_attr
        self.ref_class = ref_class
        self._selection_ref = None

        dummy_trigger = lambda item: True
        self._trigger_add_pre = triggers_add[0] or dummy_trigger
        self._trigger_add_post = triggers_add[1] or dummy_trigger
        self._trigger_del_pre = triggers_del[0] or dummy_trigger
        self._trigger_del_post = triggers_del[1] or dummy_trigger

        self.read_only = False
        self.extend(items)
        self.read_only = bool(read_only)

    @property
    def keys(self):
        """
        Return a list of the keys of items in the list.
        """
        if self.key_attr is not None:
            return [getattr(i, self.key_attr) for i in self.data]
        else:
            return list(range(len(self.data)))

    @property
    def selection_ref(self):
        """
        Return None or the original list from which this list has been
        created, as returned by :meth:`select`, :meth:`select_one` or
        :meth:`order_by`.
        """
        return self._selection_ref

    def _create_selection(self, sel_objects):
        cls = self.__class__
        selection = cls(sel_objects, key_attr=self.key_attr)
        if self._selection_ref is None:
            selection._selection_ref = self
        else:
            selection._selection_ref = self._selection_ref
        return selection

    def _get_selection_ref_or_error(self, select_one=False):
        if self._selection_ref is None:
            raise ValueError("Invalid {0}: this {0} must be "
                             "returned by the 'select' or 'select_one' method"
                             .format(self.__class__.__name__))
        elif select_one and len(self.data) != 1:
            raise ValueError("Invalid {0}: this {0} must be "
                             "returned by the 'select_one' method"
                             .format(self.__class__.__name__))
        return self._selection_ref

    def _check_read_only(self):
        if self.read_only:
            raise exceptions.NotPermittedError("read-only {0}"
                                               .format(self.__class__.__name__))

    def _check_before_add(self, item):
        if self.ref_class is None:
            self.ref_class = item.__class__

        self._check_read_only()
        if not isinstance(item, self.ref_class):
            raise ValueError(
                "cannot add a '{0}' item in a list of '{1}' items"
                .format(type(item).__name__, self.ref_class.__name__)
            )

        item_key = getattr(item, self.key_attr)
        if item_key in self.keys:
            raise ValueError(
                "a '{0}' item with key '{0}' is already in the list"
                .format(self.ref_class.__name__, item_key)
            )

    def selection_replace(self, new_item):
        """
        If this list is a selection of one object, it replace the latter
        by another object `new_item` in the original list given by
        :prop:`selection_ref`.
        """
        self._check_read_only()
        selection_ref = self._get_selection_ref_or_error(select_one=True)
        index = self.selection_index()[0]

        del selection_ref[index]
        selection_ref.insert(new_item, index=index)

    def selection_remove(self):
        """
        If this list is a selection, it removes all the corresponding items
        in the original list given by :prop:`selection_ref`.
        """
        self._check_read_only()
        selection_ref = self._get_selection_ref_or_error()
        for index in self.selection_index():
            del selection_ref[index]

    def selection_index(self):
        """
        Get the indexes of selected items in the original list given by
        :prop:`selection_ref`.
        """
        selection_ref = self._get_selection_ref_or_error()
        return [selection_ref.data.index(item) for item in self.data]

    def order_by(self, key):
        """
        Sort objects in the list.

        Parameters
        ----------
        key : string or callable
            Name of the attribute on which sorting is based, or a key function
            - which accepts any item in the set as one argument - used
            to extract a comparison key from each item of the set
            (see :func:`sorted`).

        Returns
        -------
        A new :class:`RecordList` instance.

        See Also
        --------
        :func:`sorted`
        """
        try:
            sorted_objects = sorted(self.data, key=key)
        except TypeError:
            sorted_objects = sorted(self.data,
                                    key=lambda obj: getattr(obj, key))
        return self._create_selection(sorted_objects)

    def select(self, *args, **kwargs):
        """
        Select objects in the list based on their attributes.

        Parameters
        ----------
        *args
            Any key value(s), or any callable(s) that accepts any object of
            `ref_class` as argument and that returns True or False.
        **kwargs
            Used to define simple conditions, i.e., the value of the
            attribute specified by the keyword equals the given value.

        Returns
        -------
        A new list (:class:`RecordList` instance) containing the
        selected objects.

        Notes
        -----
        If several args or kwargs are given, the returned selection satisfies
        all conditions (similar to the 'AND' operator).

        See Also
        --------
        :prop:`keys`
        :meth:`select_one`
        :meth:`select_object`
        """
        selection = self.data
        for a in args:
            if not callable(a):
                selection = filter(lambda obj: getattr(obj, self.key_attr) == a,
                                   selection)
            else:
                selection = filter(a, selection)
        for attr_name, attr_val in kwargs.items():
            selection = filter(lambda obj: getattr(obj, attr_name) == attr_val,
                               selection)
        return self._create_selection(selection)

    def select_one(self, *args, **kwargs):
        """
        Select one object based on its attributes.

        Similar to :meth:`select` but raises an error if the selection criteria
        result in other than one object.

        Returns
        -------
        A new 1-length list (:class:`RecordList` instance) containing
        the selected object.

        Raises
        ------
        :class:`pygchem.utils.exceptions.SelectionMismatchError`
            if more than one object is found or if no object is found.

        See Also
        --------
        :meth:`select`
        :meth:`select_object`
        """
        selection = self.select(*args, **kwargs)
        if len(selection.data) > 1:
            raise exceptions.SelectionMismatchError("More than one object match"
                                                    " the selection criteria")
        elif not selection.data:
            raise exceptions.SelectionMismatchError("No object match the "
                                                    "selection criteria")
        else:
            return selection

    def select_item(self, *args, **kwargs):
        """
        Select one item (object) based on its attributes.

        Similar to :meth:`select_one` but returns the object rather than a new
        :class:`RecordList`.

        See Also
        --------
        :meth:`select`
        :meth:`select_one`
        """
        selection = self.select_one(*args, **kwargs)
        obj = selection.data[0]
        return obj

    def to_list(self):
        """Return a list from this object."""
        return self.data

    def to_dict(self, ordered=False):
        """Return a (ordered) dictionary (`key_attr`: item) from this object."""
        if ordered:
            return OrderedDict(zip(self.keys, self.data))
        else:
            return dict(zip(self.keys, self.data))

    def insert(self, i, item):
        self._check_before_add(item)
        self._trigger_add_pre(item)
        self.data.insert(i, item)
        self._trigger_add_post(item)

    def append(self, item):
        self._check_before_add(item)
        self._trigger_add_pre(item)
        self.data.append(item)
        self._trigger_add_post(item)

    def extend(self, other):
        other = list(other)
        for item in other:
            self._check_before_add(item)
            self._trigger_add_pre(item)
        super(RecordList, self).extend(other)
        for item in other:
            self._trigger_add_post(item)

    def remove(self, item):
        self._check_read_only()
        self._trigger_del_pre(item)
        self.data.remove(item)
        self._trigger_del_post(item)

    def pop(self, index=-1):
        item = self.data[index]
        del self.data[index]
        return item

    def __setitem__(self, index, item):
        self._trigger_add_pre(item)
        self._check_before_add(item)
        super(RecordList, self).__setitem__(index, item)
        self._trigger_add_post(item)

    def __delitem__(self, index):
        self._check_read_only()
        item = self.data[index]
        self._trigger_del_pre(item)
        del self.data[index]
        self._trigger_del_post(item)

    def __add__(self, other):
        cls = self.__class__
        self._check_read_only()
        new_list = cls(self.data, key_attr=self.key_attr,
                       triggers_add=(self._trigger_add_pre,
                                     self._trigger_add_post),
                       triggers_del=(self._trigger_del_pre,
                                     self._trigger_del_post))
        new_list.extend(other)

    def __radd__(self, other):
        cls = self.__class__
        self._check_read_only()
        new_list = cls(other, key_attr=self.key_attr,
                       triggers_add=(self._trigger_add_pre,
                                     self._trigger_add_post),
                       triggers_del=(self._trigger_del_pre,
                                     self._trigger_del_post))
        new_list.extend(self.data)

    def __iadd__(self, other):
        other = list(other)
        for item in other:
            self._check_before_add(item)
            self._trigger_add_pre(item)
        super(RecordList, self).__iadd__(other)
        for item in other:
            self._trigger_add_post(item)

    def __mul__(self, other):
        raise exceptions.NotPermittedError('duplicating items is not allowed')

    def __imul__(self, other):
        raise exceptions.NotPermittedError('duplicating items is not allowed')

    def __str__(self):
        return "List of {0} {1}{2}:\n{3}" \
            .format(len(self), self.ref_class.__name__,
                    ' (selection)' if self._selection_ref is not None else '',
                    '\n'.join(str(obj) for obj in self.data))

    def __repr__(self):
        return "<{0}{1}: {2}>"\
            .format(self.__class__.__name__,
                    ' (selection)' if self._selection_ref is not None else '',
                    super(RecordList, self).__repr__())

    # TODO: __html__ representation


class Record(object):
    """
    Base class for record-like objects.

    Inherited by classes created with :func:`record_cls`.

    """
    __slots__ = tuple()
    _properties = tuple()
    _types = tuple()

    def __init__(self, *args):
        # convert type if needed
        for k, t, v in zip(self.__slots__, self._types, args):
            if t is not None:
                v = t(v)
            setattr(self, k, v)

    def keys(self):
        for k in self._properties:
            yield k

    def values(self):
        for k in self._properties:
            yield getattr(self, k)

    def items(self):
        for k in self._properties:
            yield (k, getattr(self, k))

    def to_dict(self, ordered=False):
        """Return a new dict which maps field names to their values."""
        kv_pairs = zip(self._properties,
                       (getattr(self, k) for k in self._properties))
        if ordered:
            return OrderedDict(kv_pairs)
        else:
            return dict(kv_pairs)

    def to_list(self):
        """Return a new list with ordered field values."""
        return list(self)

    def __len__(self):
        return len(self._properties)

    def __iter__(self):
        return self.values()

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        return setattr(self, key, value)

    def __eq__(self, other):
        compare_pairs = (getattr(self, k) == getattr(other, k)
                         for k in self._properties)
        return isinstance(other, self.__class__) and all(compare_pairs)

    def __ne__(self, other):
        return not self == other

    def __getstate__(self):
        return tuple(getattr(self, k) for k in self._properties)

    def __setstate__(self, state):
        for k, v in zip(self._properties, state):
            setattr(self, k, v)

    def __repr__(self):
        kv_pairs = zip(self._properties,
                       (getattr(self, k) for k in self._properties))
        kv_repr = ', '.join("{0}={1}".format(k, v) for k, v in kv_pairs)
        return "({0})".format(kv_repr)


def record_cls(cls_name, cls_description, fields, required_fields=(),
               verbose=False):
    """
    Create a new class designed for custom record-like objects with
    named fields.

    Parameters
    ----------
    cls_name : string
        name of the returned class.
    cls_description : string
        class short description (will appear in the class's docstring).
    fields : sequence of tuples (name, type, default, read_only, description)
        field specifications. 'name' is the name of the field, 'type' is the
        type expected for the field (should be a built-in type like bool, str,
        int, float, dict... or None if no specific type is expected), 'default'
        is the default value of the field (useless if field is required),
        read_only (bool) defines whether the field value can be modified or not
        after having created a record instance, and 'description' is the field
        description (docstring).
    required_fields : sequence of strings
        names of the fields for which value must be set explicitly.
    verbose : bool
        if True, print the code used to define the class.

    Returns
    -------
    class
        the created record-like class.

    """
    # Taken and modified from this recipe:
    # http://code.activestate.com/recipes/576555-records/ (MIT License)
    tfield = ('name', 'type', 'default', 'read_only', 'description')
    fields = tuple((dict(zip(tfield, f)) for f in fields))
    for f in fields:
        if f['type'] == str:
            f['default'] = "'{0}'".format(f['default'])
        if f['type'] is None:
            f['type'] = ''
        else:
            f['type'] = f['type'].__name__
    optional_fields = tuple(f for f in fields
                            if f['name'] not in required_fields)
    required_fields = tuple(f for f in fields if f['name'] in required_fields)
    fields = required_fields + optional_fields

    # Parse and validate the field names (prevent template injection attacks)
    field_names = tuple(map(str, (f['name'] for f in fields)))
    if not field_names:
        raise ValueError("Records must have at least one field")
    for name in (cls_name,) + field_names:
        if not min(c.isalnum() or c == '_' for c in name):
            raise ValueError("Class name and field names can only contain "
                             "alphanumeric characters and underscores: {0}"
                             .format(name))
        if iskeyword(name):
            raise ValueError("Class name and field names cannot be a keyword: "
                             "{0}".format(name))
        if name[0].isdigit():
            raise ValueError("Class name and field names cannot start with a "
                             "number: {0}".format(name))
        if name.startswith('_'):
            raise ValueError("Class name and field names cannot start with an "
                             "underscore: {0}".format(name))
    for name, count in Counter(field_names).items():
        if count > 1:
            raise ValueError("Encountered duplicate field name: {0}"
                             .format(name))

    # Create and fill-in the class template
    def indent(text, n=1):
        return os.linesep.join((n * "    ") + i for i in text.splitlines())

    parameter_template = indent(dedent("""
        {name} : {type}
            {description}
    """))
    property_template = indent(dedent("""
        @property
        def {name}(self):
            '''
            {description}.
            '''
            return self._{name}

    """))
    property_setter_template = indent(dedent("""
        @{name}.setter
        def {name}(self, value):
            self._{name} = {type}(value)

    """))

    class_template = dedent("""
        class {cls_name}(Record):
            '''
            {cls_description} (record-like).

            Parameters
            ----------{parameters_str}

            '''

            __slots__  = ({slots_str})
            _properties = {field_names}
            _types = ({types_str})

            def __init__(self, {args_kwargs_str}):
                super({cls_name}, self).__init__({init_str})

            {properties_str}
            {properties_setter_str}

            def __repr__(self):
                return '{cls_name}' + super({cls_name}, self).__repr__()

    """)

    parameters_str = ''.join(parameter_template.format(**f) for f in fields)
    slots_str = ', '.join("'_{name}'".format(**f) for f in fields)
    types_str = ', '.join(f['type'] if f['type'] else 'None' for f in fields)
    args_str = ', '.join("{name}".format(**f) for f in required_fields)
    kwargs_str = ', '.join("{name}={default}".format(**f)
                           for f in optional_fields)
    if not len(required_fields):
        args_kwargs_str = kwargs_str
    elif not len(optional_fields):
        args_kwargs_str = args_str
    else:
        args_kwargs_str = ', '.join((args_str, kwargs_str))
    init_str = ', '.join("{name}".format(**f) for f in fields)
    properties_str = ''.join(property_template.format(**f) for f in fields)
    properties_setter_str = ''.join(property_setter_template.format(**f)
                                    if not f['read_only'] else ''
                                    for f in fields)
    class_str = class_template.format(**locals())

    # Execute the class string in a temporary namespace
    namespace = {}
    try:
        exec "from pygchem.utils.data_structures import Record" in namespace
        exec class_str in namespace
        if verbose:
            print(class_str)
    except SyntaxError, e:
        raise SyntaxError(e.message + ':\n' + class_str)
    cls = namespace[cls_name]
    #cls.__init__.im_func.func_defaults = init_defaults

    # For pickling to work, the __module__ variable needs to be set to the frame
    # where the named tuple is created.  Bypass this step in environments where
    # sys._getframe is not defined (Jython for example).
    if hasattr(sys, '_getframe') and sys.platform != 'cli':
        cls.__module__ = sys._getframe(1).f_globals['__name__']

    return cls

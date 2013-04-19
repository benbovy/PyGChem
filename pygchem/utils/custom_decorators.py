# -*- coding: utf-8 -*-

# module pygchem.utils.custom_decorators
# pygchem: Python interface for GEOS-Chem Chemistry Transport Model
#
# Copyright (C) 2012 Benoit Bovy
# see license.txt for more details
#
# Last modification: 11/2012

"""
Module that defines customized decorators

Authors: Gribouillis
"""

from functools import partial

class mixedmethod(object):
    """
    This decorator mutates a function defined in a class into a 'mixed' class 
    and instance method.
    
    Example
    -------
    class Spam:

    @mixedmethod
    def egg(self, cls, *args, **kwargs):
        if self is None:
            pass # executed if egg was called as a class method
        else:
            pass # executed if egg was called as an instance method

    The decorated methods need 2 implicit arguments: self and cls, the former 
    being None when there is no instance in the call. 
    This follows the same rule as __get__ methods in python's 
    descriptor protocol.
    """
    def __init__(self, func):
        self.func = func
    def __get__(self, instance, cls):
        return partial(self.func, instance, cls)

class classproperty(property):
    """
    A class property
    
    class Foo(object):
        _var = 5

        @ClassProperty
        @classmethod
        def var(cls):
            return cls._var

        @var.setter
        @classmethod
        def var(cls, value):
            cls._var = value
    """
    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()

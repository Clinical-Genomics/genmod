#!/usr/bin/env python
# encoding: utf-8
"""
is_number.py

Created by Måns Magnusson on 2013-04-09.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""
from __future__ import print_function


def is_number(number):
    """
    Returns true if the input is a number or False otherwise
    
    Arguments:
        number (obj): The object that should be checked
    
    """
    try:
        float(number)
        return True
    except ValueError:
        pass
    return False


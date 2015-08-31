#!/usr/bin/env python
# encoding: utf-8
"""
pair_generator.py

Class that takes a list of objects and return all unordered pairs as a generator.

If only one object? Raise Exception
 
Created by Måns Magnusson on 2013-03-01.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function

import collections

def generate_pairs(objects):
    """
    Yields all unordered pairs as tuples from the list of objects
    
    Arguments:
        list_of_objects (iterator):
    """
    if not isinstance(objects, collections.Iterable):
        raise SyntaxError("objects has to be iterable. objects: {0}".format(
            objects
        ))
    if len(objects) < 2:
        #TODO raise a proper exception here
        raise SyntaxError('List must include at least 2 objects!."\
                        " objects: {0}'.format(objects))
        
    for i in range(len(objects)-1):
        for j in range(i+1, len(objects)):
            yield (objects[i], objects[j])
    


#!/usr/bin/env python
# encoding: utf-8
"""
pair_generator.py

Class that takes a list of objects and return all unordered pairs as a generator.

If only one object? Raise Exception
 
Created by Måns Magnusson on 2013-03-01.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse

class Pair_Generator(object):
    """Yeilds all unordered pairs from a list of objects as tuples, like (obj_1, obj_2)"""
    def __init__(self, list_of_objects):
        super(Pair_Generator, self).__init__()
        if len(list_of_objects) < 2:
            #TODO raise a proper exception here
            raise SyntaxError('List must include at least 2 objects!. List: %s' % list_of_objects)
        self.list_of_objects = list_of_objects
    
    def generate_pairs(self):
        """Yields all unordered pairs as tuples from the list of objects"""
        for i in range(len(self.list_of_objects)-1):
            for j in range(i+1, len(self.list_of_objects)):
                yield (self.list_of_objects[i], self.list_of_objects[j])
    

def main():
    parser = argparse.ArgumentParser(description="Generate all unordered pairs of a sequence.")
    parser.add_argument('-seq','--sequence', nargs='+', default=['a', 'b', 'c', 'd'], help='A sequence of objects.')
    args = parser.parse_args()
    
    for pairs in Pair_Generator(args.sequence).generate_pairs():
        print(pairs)


if __name__ == '__main__':
    main()


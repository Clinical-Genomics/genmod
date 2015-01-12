#!/usr/bin/env python
# encoding: utf-8
"""
warning.py

Prints a warning message to stderr

Created by Måns Magnusson on 2014-08-22.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

import sys
import os

def warning(*objs):
    """Prints the warning messages to std err"""
    print("WARNING: ", *objs, file=sys.stderr)
    
def main():
    print('This is a warning!')

if __name__ == '__main__':
    main()


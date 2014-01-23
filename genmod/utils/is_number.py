#!/usr/bin/env python
# encoding: utf-8
"""
is_number.py

Created by Måns Magnusson on 2013-04-09.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os

def is_number(number):
    """Returns true if the string is a number or False otherwise"""
    if type(number) == type(1) or type(number) == type(0.1) or type(number) == type(''):
        try:
            float(number)
            return True
        except ValueError:
            return False
    else:
        return False

def main():
    number = sys.argv[1]
    print is_number(number)


if __name__ == '__main__':
    main()


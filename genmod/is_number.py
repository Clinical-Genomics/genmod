#!/usr/bin/env python
# encoding: utf-8
"""
is_number.py

Created by Måns Magnusson on 2013-04-09.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse

def is_number(number):
    """Returns true if the string is a number or False otherwise"""
    if type(number) == type(1) or type(number) == type(0.1) or type(number) == type('') or type(u''):
        try:
            float(number)
            return True
        except ValueError:
            return False
        except TypeError:
            return False
    else:
        return False

def main():
    parser = argparse.ArgumentParser(description="Annotate genetic models in variant files..")
    
    parser.add_argument('input', nargs=1, 
        help='give an input to see if it is a number.'
    )
    args = parser.parse_args()
    
    number = args.input[0]
    
    print(is_number(number))


if __name__ == '__main__':
    main()


#!/usr/bin/env python
# encoding: utf-8
"""
test_is_number.py

Test the is_number function. Should return true if a string is a number false otherwise.

Created by Måns Magnusson on 2013-05-31.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from genmod.utils.is_number import is_number

def test_integer():
	"""Try if is_number behaves correct when given an integer."""
	my_integer_string = '1'
	my_integer = 1
	assert is_number(my_integer_string)
	assert is_number(my_integer)

def test_float():
	"""Try is is_number behaves correct when given a float."""
	my_float = 0.22
	my_float_string = '0.22'
	assert is_number(my_float_string)
	assert is_number(my_float)

def test_zero():
	"""Try how is_number behaves when given zero."""
	my_zero = 0
	my_zero_string = '0'
	assert is_number(my_zero)
	assert is_number(my_zero_string)

def test_negative_number():
	"""Try how is_number behaves when given a negative number."""
	my_negative = -14
	my_negative_string = '-12'
	assert is_number(my_negative)
	assert is_number(my_negative_string)

def test_non_number():
	"""Try if is_number behaves correct when not given a number."""
	my_empty_string = ''
	my_minus = '-'
	my_string = 'g'
	my_none = None
	assert not is_number(my_none)
	assert not is_number(my_empty_string)
	assert not is_number(my_string)
	assert not is_number(my_minus)

def main():
	pass


if __name__ == '__main__':
	main()


#!/usr/bin/env python
# encoding: utf-8
"""
test_pair_generator.py

Test for the variant class.

Created by Måns Magnusson on 2013-02-28.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import pytest
from genmod.utils import pair_generator

class TestPairGenerator(object):
    """Test class for testing how the Pair_Generator object behave"""
    
    def setup_class(self):
        """Setup a simple family with family id 1, sick daughter id 1, healthy father id 2, healthy mother id 3"""
        list_of_objects = ['1', '2', '3', 'f']
        self.my_pair_generator = pair_generator.Pair_Generator(list_of_objects)
    
    def test_generator(self):
        """Test if the pairs are generated in a correct way"""
        my_pairs = self.my_pair_generator.generate_pairs()
        assert next(my_pairs) == ('1', '2')
        assert next(my_pairs) == ('1', '3')
        assert next(my_pairs) == ('1', 'f')
        assert next(my_pairs) == ('2', '3')
        assert next(my_pairs) == ('2', 'f')
        assert next(my_pairs) == ('3', 'f')
        with pytest.raises(StopIteration):
            next(my_pairs)
    

def main():
    pass


if __name__ == '__main__':
    main()
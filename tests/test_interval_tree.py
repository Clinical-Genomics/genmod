#!/usr/bin/env python
# encoding: utf-8
"""
test_interval_tree.py

Test so that the interval tree behave as expected

Created by Måns Magnusson on 2014-03-12.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import tempfile
import pytest

from genmod.utils.interval_tree import IntervalTree


class TestAnnotationParser(object):
    """Test class for testing how the genetic models behave with combinations of compound variants."""

    def setup_class(self):
        """Setup different interval trees and check if they behave correct"""
        # Setup family with sick kid, sick father and healthy mother:
        # This should behave such that pos 1 and 2 are covered by the interval:
        smallest_interval = [1,1,'id_01']
        self.simplest_tree = IntervalTree([smallest_interval],1, 2)
        small_interval = [1,2,'id_02']
        self.simple_tree = IntervalTree([smallest_interval, small_interval],1, 2)
        interval_3 = [10,20,'id_03']
        interval_4 = [16,40,'id_04']
        self.tree = IntervalTree([interval_3, interval_4],1, 50)
         
    
    def test_basic_tree_functions(self):
        """Check so that the tree class behave as expected."""
        with pytest.raises(SyntaxError):
            tree = IntervalTree([[2,'id01']],1,2)
        with pytest.raises(SyntaxError):
            self.simplest_tree.find_range([1])
        with pytest.raises(ValueError):
            self.simplest_tree.find_range([1, 'a'])
        
    def test_simplest_tree(self):
        """Test interval_tree with one position."""
        assert self.simplest_tree.find_range([1,1]) == ['id_01']
        assert self.simplest_tree.find_range([2,2]) == []
        
    def test_simple_tree(self):
        """Test interval_tree with one position."""
        assert set(self.simple_tree.find_range([1,1])) == set(['id_01', 'id_02'])
        assert set(self.simple_tree.find_range([1,2])) == set(['id_01', 'id_02'])
        assert set(self.simple_tree.find_range([2,2])) == set(['id_02'])
        assert set(self.simple_tree.find_range([3,3])) == set([])

    def test_tree(self):
        """Test interval_tree with one position."""
        assert set(self.tree.find_range([1,1])) == set([])
        assert set(self.tree.find_range([10,12])) == set(['id_03'])
        assert set(self.tree.find_range([15,15])) == set(['id_03'])
        assert set(self.tree.find_range([16,16])) == set(['id_03', 'id_04'])
        assert set(self.tree.find_range([12,22])) == set(['id_03', 'id_04'])
        assert set(self.tree.find_range([30,30])) == set(['id_04'])
        assert set(self.tree.find_range([30,30])) == set(['id_04'])
        assert set(self.tree.find_range([42,42])) == set([])

def main():
    pass


if __name__ == '__main__':
    main()


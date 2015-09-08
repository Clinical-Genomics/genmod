#!/usr/bin/env python
# encoding: utf-8
"""
score_function.py

Create a score function

Created by Måns Magnusson on 2015-09-08.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import logging

from intervaltree import IntervalTree

class ScoreFunction(object):
    """Class for holding score functions"""
    def __init__(self, match_type, string_dict={}):
        super(ScoreFunction, self).__init__()
        self.match_type = match_type #['integer','float','flag','character','string']
        self.string_dict = string_dict
        self.interval_tree = IntervalTree()
        self.not_reported_score = 0
        self.reported_score = 0 # only for 'flag'
        
    def add_interval(self, lower, upper, score):
        """Add an interval to the score function
        
            Args:
                lower (int,float): The lower bound of the interval
                upper (int,float): The upper bound of the interval
                score (int,float): The score of the interval
        """
        self.interval_tree[lower:upper] = score
        
        
    def get_score(self, value):
        """Take a value and return a score
            
            If value is None we return the not_reported score
            If value is not None but does not have a rule we return 0
            If Score function is a string comparison we match the string
            If value is a number (float or int):
                if operator is equal we return the number
                else return data of interval
        """
        score = 0
        if not value:
            score = self.not_reported_score
        elif self.match_type == 'string':
            score = self.string_dict.get(value, 0)
        elif self.match_type == 'flag':
            # Here we have a value so we give the reported score
            score = self.reported_score
        elif self.match_type in ['integer', 'float']:
            for interval in self.interval_tree[value]:
                score = interval.data
        
        return score

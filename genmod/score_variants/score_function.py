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
    def __init__(self, match_type, equal=False):
        super(ScoreFunction, self).__init__()
        self.logger = logging.getLogger(__name__)
        self.logger.debug("Initializing match_type to:{0}".format(match_type))
        self.match_type = match_type #['integer','float','flag','character','string']
        self.logger.debug("Initializing string_dict to:{}")
        self._string_dict = {}
        self.logger.debug("Initializing interval_tree")
        self._interval_tree = IntervalTree()
        self.logger.debug("Initializing value_dict")
        self._value_dict = {}
        self.logger.debug("Initializing not_reported_score to 0")
        self._not_reported_score = 0
        self.logger.debug("Initializing reported_score to 0")
        self._reported_score = 0 # only for 'flag'
        # If the score is the same as the value found:
        self.logger.debug("Initializing equal to {0}".format(equal))
        self._equal = equal
        
    def add_interval(self, lower, upper, score):
        """Add an interval to the score function
        
            Args:
                lower (int,float): The lower bound of the interval
                upper (int,float): The upper bound of the interval
                score (int,float): The score of the interval
        """
        self.logger.debug("Adding interval {0} to score function".format(
            ','.join([str(lower), str(upper), str(score)])
        ))
        ##TODO Check if intervals overlap
        self._interval_tree[lower:upper] = score
        
        return
    
    def add_string_rule(self, key, score):
        """Add the score for a string match
        
            Args:
                key (str): The string that should be matched
                score (int): The score for the match
            
        """
        self.logger.debug("Adding string {0} with score {1} to string_dict".format(
            key, str(score))
        )
        self._string_dict[key.lower()] = score
        return

    def add_value(self, value, score):
        """Add the score for a value match
        
            Args:
                value (number): The number that should be matched
                score (int): The score for the match
            
        """
        self.logger.debug("Adding value {0} with score {1} to value_dict".format(
            value, str(score))
        )
        self._value_dict[str(value)] = score
        return
        
    def get_score(self, value):
        """Take a value and return a score
            
           - If value is None we return the not_reported score
           - If value is NOT None but does not have a rule we return 0
           - If Score function is a string comparison we match the string
           - If value is a number (float or int):
                - if operator is equal we return the number
                - else return data of interval
            
            Args:
                value (str): The value that we want to find the score for
            
            Return:
                score (number): The score for this value
        """
        score = 0
        
        if not value:
            self.logger.debug("No value found set score to not reported score")
            score = self._not_reported_score
        
        # Here we know there is a value 
        elif self.match_type == 'flag':
            self.logger.debug("Flag found set score reported score")
            score = self._reported_score
        
        elif self.match_type in ['string', 'char']:
            score = self._string_dict.get(value.lower(), 0)
        
        # We know that match type must be any of integer or float
        else:
            if self.match_type == 'float':
                try:
                    value = float(value)
                except ValueError:
                    raise ValueError("Value has to be a number")
            else:
                try:
                    value = int(value)
                except ValueError:
                    raise ValueError("Value has to be a number")
            
            if self._equal:
                score = value
            
            else:
                if self._value_dict:
                    score = float(self._value_dict.get(str(value), 0))
                    self.logger.debug("Got score from value dict")
                else:
                    #There should only be one interval matching
                    ##TODO Check if intervals overlap
                    for interval in self._interval_tree[value]:
                        score = interval.data
                        self.logger.debug("Got score from interval tree")
        
        # For now we only allow integers as score
        ## TODO fix this ugly solution
        try:
            score = int(score)
        except error as e:
            score = int(float(score))
            
        
        return score
    
    def set_not_reported(self, value):
        """Set the not reported score
        
        Args:
            value (int, float): The not reported score
        
        """
        self.logger.debug("Setting not_reported_score to {0}".format(value))
        self._not_reported_score = float(value)
        return

    def set_reported(self, value):
        """Set the reported score
        
        Args:
            value (int, float): The reported score
        """
        self.logger.debug("Setting reported_score to {0}".format(value))
        self._reported_score = float(value)
        return
    
    def set_equal(self):
        """Set _equal to True
        """
        self.logger.debug("Setting equal to True")
        self._equal = True
        return

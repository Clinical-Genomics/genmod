#!/usr/bin/env python
# encoding: utf-8
"""
test_annotation_parser.py

Test so that the annotation parser behaves as expected with the different file formats

Created by Måns Magnusson on 2014-03-12.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import tempfile

from genmod.utils.annotation_parser import AnnotationParser


class TestAnnotationParser(object):
    """Test class for testing how the genetic models behave with combinations of compound variants."""

    def setup_class(self):
        """Setup a simple family with family id 1, sick son id 1,
         healthy father id 2, healthy mother id 3"""
        # Setup family with sick kid, sick father and healthy mother:
        bed_header = '#\t'
        #TODO write proper test for these
        
    def test_bed(self):
        """Test the bed file parser."""
        #TODO write bed tests
    
    def test_ccds(self):
        """Test the ccds file parser."""
        #TODO write ccds tests
        
    def test_gtf(self):
        """Test gtf file parser."""
        #TODO write bed tests
        
    def test_ref_gene(self):
        """Test ref_gene parser."""
        #TODO write bed tests
    

def main():
    pass


if __name__ == '__main__':
    main()


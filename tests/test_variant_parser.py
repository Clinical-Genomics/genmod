#!/usr/bin/env python
# encoding: utf-8
"""
test_variant_parser.py

Test the so that the variant parser behave as suspected.


Created by Måns Magnusson on 2014-03-04.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

import sys
import os

from pprint import pprint as pp

from ped_parser import family, individual
from genmod.vcf import vcf_parser
from genmod.vcf import vcf_header


class TestVariantParser(object):
    """Test how how the variant parser behaves."""

    def setup_class(self):
        """Setup a vcf file with some variants and a interval tree with features."""
        pass
        #TODO make tests        
        

def main():
    pass


if __name__ == '__main__':
    main()


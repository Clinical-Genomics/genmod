from __future__ import print_function

import sys
import re

from logging import getLogger

if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict


class HeaderParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self):
        super(HeaderParser, self).__init__()
        self.logger = getLogger(__name__)
        self.info_lines=[]
        self.info_dict=OrderedDict()
        #This is a dictionary cantaining specific information about the info fields
        #It will have info name as key and then another dictionary with ID, Number, Type and Description
        self.extra_info = {}
        
        self.filter_lines=[]
        self.filter_dict=OrderedDict()
        
        self.contig_lines=[]
        self.contig_dict=OrderedDict()
        
        self.format_lines=[]
        self.format_dict=OrderedDict()
        
        self.alt_lines=[]
        self.alt_dict=OrderedDict()
        
        self.other_lines=[]
        self.other_dict=OrderedDict()
        
        self.header=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
        self.header_keys={'info' : ['ID', 'Number', 'Type', 'Description'], 
                            'form' : ['ID', 'Number', 'Type', 'Description'], 
                            'filt' : ['ID', 'Description'],
                            'alt' : ['ID', 'Description'],
                            'contig' : ['ID', 'length']}
        self.fileformat = None
        self.filedate = None
        self.reference = None
        self.phasing = None
        self.source = None
        self.line_counter = 0
        self.individuals = []
        self.vep_columns = []
        self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AGR]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.contig_pattern = re.compile(r'''\#\#contig=<
            ID=(?P<id>[^,]+),
            .*
            length=(?P<length>-?\d+)
            .*
            >''', re.VERBOSE)
        self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AGR]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.alt_pattern = re.compile(r'''\#\#ALT=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')
    
    def parse_meta_data(self, line):
        """Parse a vcf metadataline"""
        line = line.rstrip()
        line_info = line[2:].split('=')
        match = False
        
        if line_info[0] == 'fileformat':
            try:
                self.fileformat = line_info[1]
            except IndexError:
                raise SyntaxError("fileformat must have a value")
        
        elif line_info[0] == 'INFO':
            match = self.info_pattern.match(line)
            if not match:
                raise SyntaxError("One of the INFO lines is malformed:{0}".format(line))
            
            matches = [
                match.group('id'), match.group('number'), 
                match.group('type'), match.group('desc')
            ]
            
            # extra_info is a dictionary to check the metadata about the INFO values:
            self.extra_info[matches[0]] = dict(
                zip(self.header_keys['info'][1:], matches[1:])
            )
            
            info_line = dict(list(zip(self.header_keys['info'],matches)))
            
            if len(info_line['Description'].split('Format:')) > 1:
                info_line['Format'] = [
                    info.strip() for info in info_line['Description'].split('Format:')
                ][-1]
            self.info_lines.append(info_line)
            
            # Store the vep columns:
            if info_line['ID'] == 'CSQ':
                self.vep_columns = info_line.get('Format', '').split('|')
            
            self.info_dict[match.group('id')] = line
        
        elif line_info[0] == 'FILTER':
            match = self.filter_pattern.match(line)
            if not match:
                raise SyntaxError("One of the FILTER lines is malformed: {0}".format(line))
            matches = [match.group('id'), match.group('desc')]
            self.filter_lines.append(dict(
                list(zip(self.header_keys['filt'],matches)))
            )
            self.filter_dict[match.group('id')] = line
        
        elif line_info[0] == 'contig':
            match = self.contig_pattern.match(line)
            if not match:
                print()
                raise SyntaxError("One of the contig lines is malformed: {0}".format(line))
            
            matches = [match.group('id'), match.group('length')]
            self.contig_lines.append(dict(
                list(zip(self.header_keys['contig'],matches)))
            )
            self.contig_dict[match.group('id')] = line
        
        elif line_info[0] == 'FORMAT':
            match = self.format_pattern.match(line)
            if not match:
                raise SyntaxError("One of the FORMAT lines is malformed: {0}".format(line))
            
            matches = [
                match.group('id'), match.group('number'), 
                match.group('type'), match.group('desc')
            ]
            self.format_lines.append(dict(
                list(zip(self.header_keys['form'],matches)))
            )
            self.format_dict[match.group('id')] = line
        
        elif line_info[0] == 'ALT':
            match = self.alt_pattern.match(line)
            if not match:
                raise SyntaxError("One of the ALT lines is malformed: {0}".format(line))
            
            matches = [match.group('id'), match.group('desc')]
            self.alt_lines.append(dict(
                list(zip(self.header_keys['alt'],matches)))
            )
            self.alt_dict[match.group('id')] = line
        
        else:
            match = self.meta_pattern.match(line)
            if not match:
                raise SyntaxError("One of the meta data lines is malformed: {0}".format(line))
            
            self.other_lines.append({match.group('key'): match.group('val')})
            self.other_dict[match.group('key')] = line
    
    def parse_header_line(self, line):
        """docstring for parse_header_line"""
        self.header = line[1:].rstrip().split('\t')
        if len(self.header) < 9:
            self.header = line[1:].rstrip().split()
        self.individuals = self.header[9:]
    
    def print_header(self):
        """Returns a list with the header lines if proper format"""
        lines_to_print = []
        lines_to_print.append('##fileformat='+self.fileformat)
        if self.filedate:
            lines_to_print.append('##fileformat='+self.fileformat)
            
        for filt in self.filter_dict:
            lines_to_print.append(self.filter_dict[filt])
        for form in self.format_dict:
            lines_to_print.append(self.format_dict[form])
        for info in self.info_dict:
            lines_to_print.append(self.info_dict[info])
        for contig in self.contig_dict:
            lines_to_print.append(self.contig_dict[contig])
        for alt in self.alt_dict:
            lines_to_print.append(self.alt_dict[alt])
        for other in self.other_dict:
            lines_to_print.append(self.other_dict[other])
        lines_to_print.append('#'+ '\t'.join(self.header))
        return lines_to_print
    

    def add_fileformat(self, fileformat):
        """
        Add fileformat line to the header.
        
        Arguments:
            fileformat (str): The id of the info line
        
        """
        self.fileformat = fileformat
        self.logger.info("Adding fileformat to vcf: {0}".format(fileformat))
        return

    def add_meta_line(self, key, value):
        """
        Adds an arbitrary metadata line to the header.
        
        This must be a key value pair
        
        Arguments:
            key (str): The key of the metadata line
            value (str): The value of the metadata line
        
        """
        meta_line = '##{0}={1}'.format(
            key, value
        )
        self.logger.info("Adding meta line to vcf: {0}".format(meta_line))
        self.parse_meta_data(meta_line)
        return

    def add_info(self, info_id, number, entry_type, description):
        """
        Add an info line to the header.
        
        Arguments:
            info_id (str): The id of the info line
            number (str): Integer or any of [A,R,G,.]
            entry_type (str): Any of [Integer,Float,Flag,Character,String]
            description (str): A description of the info line
        
        """
        info_line = '##INFO=<ID={0},Number={1},Type={2},Description="{3}">'.format(
            info_id, number, entry_type, description
        )
        self.logger.info("Adding info line to vcf: {0}".format(info_line))
        self.parse_meta_data(info_line)
        return

    def add_filter(self, filter_id, description):
        """
        Add a filter line to the header.
        
        Arguments:
            filter_id (str): The id of the filter line
            description (str): A description of the info line
        
        """
        filter_line = '##FILTER=<ID={0},Description="{1}">'.format(
            filter_id, description
        )
        self.logger.info("Adding filter line to vcf: {0}".format(filter_line))
        self.parse_meta_data(filter_line)
        return

    def add_format(self, format_id, number, entry_type, description):
        """
        Add a format line to the header.
        
        Arguments:
            format_id (str): The id of the format line
            number (str): Integer or any of [A,R,G,.]
            entry_type (str): Any of [Integer,Float,Flag,Character,String]
            description (str): A description of the info line
        
        """
        format_line = '##FORMAT=<ID={0},Number={1},Type={2},Description="{3}">'.format(
            format_id, number, entry_type, description
        )
        self.logger.info("Adding format line to vcf: {0}".format(format_line))
        self.parse_meta_data(format_line)
        return

    def add_alt(self, alt_id, description):
        """
        Add a alternative allele format field line to the header.
        
        Arguments:
            alt_id (str): The id of the alternative line
            description (str): A description of the info line
        
        """
        alt_line = '##ALT=<ID={0},Description="{1}">'.format(
            alt_id, description
        )
        self.logger.info("Adding alternative allele line to vcf: {0}".format(alt_line))
        self.parse_meta_data(alt_line)
        return

    def add_contig(self, contig_id, length):
        """
        Add a contig line to the header.
        
        Arguments:
            contig_id (str): The id of the alternative line
            length (str): A description of the info line
        
        """
        contig_line = '##contig=<ID={0},length={1}>'.format(
            contig_id, length
        )
        self.logger.info("Adding contig line to vcf: {0}".format(contig_line))
        self.parse_meta_data(contig_line)
        return

    def add_contig(self, contig_id, length):
        """
        Add a contig line to the header.
        
        Arguments:
            contig_id (str): The id of the alternative line
            length (str): A description of the info line
        
        """
        contig_line = '##contig=<ID={0},length={1}>'.format(
            contig_id, length
        )
        self.logger.info("Adding contig line to vcf: {0}".format(contig_line))
        self.parse_meta_data(contig_line)
        return

    def add_version_tracking(self, info_id, version, date, command_line=''):
        """
        Add a line with information about which software that was run and when 
        to the header.
        
        Arguments:
            info_id (str): The id of the info line
            version (str): The version of the software used
            date (str): Date when software was run
            command_line (str): The command line that was used for run
        
        """
        other_line = '##Software=<ID={0},Version={1},Date="{2}",CommandLineOptions="{3}">'.format(
            info_id, version, date, command_line) 
        self.other_dict[info_id] = other_line
        return

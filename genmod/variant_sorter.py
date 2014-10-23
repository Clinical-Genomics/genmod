#!/usr/bin/env python
# encoding: utf-8
"""
file_sorter.py

Sort a variant file based on position or rank score.


Created by Tomasz Beiruta on 2005-05-31.

Modified by Måns Magnusson on 2014-01-14.

"""

from __future__ import print_function, unicode_literals

import sys
import os
import click
import re
from tempfile import NamedTemporaryFile
from codecs import open

from genmod.is_number import is_number


def get_cadd_score(variant_line):
    """Return the cadd score if any, as a float."""
    cadd_score = 0.0
    for info in variant_line.split('\t')[7].split(';'):
        info = info.split('=')
        # print(info)
        if info[0] == 'CADD':
            cadd_score = max([float(cscore) for cscore in info[-1].split(',')])
    # print(cadd_score)
    return cadd_score


class FileSort(object):
    def __init__(self, inFile, mode='pos', outfile=None, splitSize=20, silent=False):
        """ split size (in MB) """
        self._inFile = inFile
        self._silent = silent
        self._mode = mode
                
        self._outFile = outfile
        self._splitSize = splitSize * 1000000
                
        if mode == 'cadd':
            # this will only work if variants have cadd scores so it needs to be fixed before
            # self._getKey = lambda variant_line: (float(re.findall(r"CADD=(\d+.\d+)", variant_line)[0]))
            self._getKey = lambda variant_line: get_cadd_score(variant_line)
        else:
            self._getKey = lambda variant_line: (int(variant_line.split('\t')[1]))
        
    
    def sort(self):
        
        files = self._splitFile()

        if files is None:
            """ file size <= self._splitSize """            
            self._sortFile(self._inFile, outFile=self._outFile, ready_to_print=True)
            return

        for fn in files:
            self._sortFile(fn)
            
        self._mergeFiles(files)
        self._deleteFiles(files)

        
    def _sortFile(self, fileName, outFile=None, ready_to_print=False):
        try:
            lines = open(fileName, mode='r', encoding='utf-8').readlines()
        except TypeError:
            lines = open(fileName.name, mode='r', encoding='utf-8').readlines()
        get_key = self._getKey
        data = [(get_key(line), line) for line in lines if line!='']
        if self._mode == 'cadd':
            data.sort(reverse=True)
        else:
            data.sort()
        lines = [line[1] for line in data]
        if ready_to_print:
            if outFile:
                try:
                    open(outFile, mode='a', encoding='utf-8').write(''.join(lines))
                except TypeError:
                #If we deal with temporary files:
                    open(outFile.name, mode='a', encoding='utf-8').write(''.join(lines))
                
            else:
                if not self._silent:
                    for line in lines:
                        print(line.rstrip())
        else:
        # In this case the temporary files are over witten.
            try:
                f = open(fileName, mode='w', encoding='utf-8')
            except TypeError:
            #If we deal with temporary files:
                f = open(fileName.name, mode='w', encoding='utf-8')
            f.write(''.join(lines))
    
    def _splitFile(self):
        totalSize = os.path.getsize(self._inFile)
        if totalSize <= self._splitSize:
            # do not split file, the file isn't so big.
            return None

        fileNames = []
        with open(self._inFile, mode='r', encoding='utf-8') as f:
            size = 0
            lines = []
            for line in f:
                size += len(line)
                lines.append(line)
                if size >= self._splitSize:
                    temp_file = NamedTemporaryFile(delete=False)
                    temp_file.close()
                    with open(temp_file.name, mode='w', encoding='utf-8') as f:
                        fileNames.append(temp_file.name)
                        f.write(''.join(lines))
                        del lines[:]
                        size = 0
                                                       
            if size > 0:
                temp_file = NamedTemporaryFile(delete=False)
                temp_file.close()
                with open(temp_file.name, mode='w', encoding='utf-8') as f:
                    fileNames.append(temp_file.name)
                    f.write(''.join(lines))
            
        return fileNames

    def _mergeFiles(self, files):
        try:
            files = [open(f, mode='r', encoding='utf-8') for f in files]
        except TypeError:
            files = [open(f.name, mode='r', encoding='utf-8') for f in files]
            
        lines = []
        keys = []
        
        for f in files:
            l = f.readline()  
            lines.append(l)
            keys.append(self._getKey(l))
        
        buff = []
        buffSize = self._splitSize/2
        append = buff.append
        if self._outFile:
            output = open(self._outFile, mode='a', encoding = 'utf-8')
        try:
            key = max(keys)
            index = keys.index(key)
            get_key = self._getKey
            while 1:
                while key == max(keys):
                    append(lines[index])
                    if len(buff) > buffSize:
                        if self._outFile:
                            output.write(''.join(buff))
                        else:
                            if not self._silent:
                                print(''.join(buff))
                        del buff[:]
                            
                    line = files[index].readline()
                    if not line:
                        files[index].close()
                        del files[index]
                        del keys[index]
                        del lines[index]
                        break
                    key = get_key(line)
                    keys[index] = key
                    lines[index] = line
        
                if len(files)==0:
                    break
                # key != min(keys), see for new index (file)
                key = max(keys)
                index = keys.index(key)

            if len(buff)>0:
                if self._outFile:
                    output.write(''.join(buff))
                else:
                    if not self._silent:
                        print(''.join(buff))
        finally:
            if self._outFile:
                output.close()
    
    def _deleteFiles(self, files):   
        try:
            for fn in files:
                os.remove(fn)
        except TypeError:
            for fn in files:
                os.remove(fn.name)
    

@click.command()
@click.argument('infile', 
                    nargs=1, 
                    type=click.Path(exists=True),
                    metavar='<vcf_file> or "-"'
)
@click.option('-cadd' ,'--cadd', 
                    is_flag=True,
                    help='Sort the variants on their cadd values.'
)
@click.option('-o' ,'--out', 
                    type=click.Path(exists=False),
                    help='Specify the path to the outfile.'
)

def cli(infile, out, cadd):
    temp_file = NamedTemporaryFile(delete=False)
    temp_file.close()
    with open(infile, mode='r', encoding = 'utf-8') as f:
        with open(temp_file.name, mode='w', encoding = 'utf-8') as g:
            for line in f:
                if not line.startswith('#'):
                    g.write(line)
    print('no errors')
    fs = FileSort(temp_file.name, args.outfile[0])
    fs.sort()


if __name__ == '__main__':
    cli()

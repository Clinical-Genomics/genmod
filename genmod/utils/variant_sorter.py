#!/usr/bin/env python
# encoding: utf-8
"""
file_sorter.py

Sort a variant file based on position or rank score.


Created by Tomasz Beiruta on 2005-05-31.
Modified by Måns Magnusson on 2014-01-14.

"""

import sys
import os
import argparse
from tempfile import NamedTemporaryFile

from genmod.utils.is_number import is_number


class FileSort(object):
    def __init__(self, inFile, outFile=None, splitSize=20, silent=False):
        """ split size (in MB) """
        self._inFile = inFile
        self._silent = silent
                
        self._outFile = outFile
                    
        self._splitSize = splitSize * 1000000
                
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
            lines = open(fileName).readlines()
        except TypeError:
            lines = open(fileName.name).readlines()
        get_key = self._getKey
        data = [(get_key(line), line) for line in lines if line!='']
        data.sort()
        lines = [line[1] for line in data]
        if ready_to_print:
            if outFile:
                open(outFile, 'a').write(''.join(lines))
            else:
                if not self._silent:
                    print ''.join(lines)
        else:
        # In this case the temporary files are over witten.
            try: 
                f = open(fileName, 'w')
            except TypeError:
                f = open(fileName.name, 'w')
            f.write(''.join(lines))
    
    def _splitFile(self):
        totalSize = os.path.getsize(self._inFile)
        if totalSize <= self._splitSize:
            # do not split file, the file isn't so big.
            return None

        fileNames = []            
        with open(self._inFile, 'r+b') as f:
            size = 0
            lines = []
            for line in f:
                if not is_number(line.rstrip().split('\t')[1]):
                    print 'hej',line
                size += len(line)
                lines.append(line)
                if size >= self._splitSize:
                    tmpFile = NamedTemporaryFile(delete=False)
                    fileNames.append(tmpFile.name)
                    tmpFile.write(''.join(lines))
                    tmpFile.close()
                    del lines[:]
                    size = 0
                                                       
            if size > 0:
                tmpFile = NamedTemporaryFile(delete=False)
                fileNames.append(tmpFile.name)
                tmpFile.write(''.join(lines))
                tmpFile.close()
            return fileNames

    def _mergeFiles(self, files):
        try:
            files = [open(f, 'r+b') for f in files]
        except TypeError:
            files = [open(f.name, 'r+b') for f in files]
            
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
            output = open(self._outFile,'a')
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
                                print ''.join(buff)
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
                        print ''.join(buff)
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
    


def main():
    parser = argparse.ArgumentParser(description="Check files.")
    parser.add_argument('infile', type=str, nargs=1, help='Specify the path to the file of interest.')
    parser.add_argument('-out', '--outfile', type=str, nargs=1, default=[None], help='Specify the path to the outfile.')
    args = parser.parse_args()
    infile = args.infile[0]
    new_file = NamedTemporaryFile(delete=False)
    with open(infile, 'rb') as f:
        for line in f:
            if not line.startswith('#'):
                new_file.write(line)
    for line in new_file.readlines():
        if not is_number(line.rstrip().split('\t')[-1]):
            print 'du', line
    print 'no errors'
    fs = FileSort(new_file, args.outfile[0])
    fs.sort()
                    
                
             


if __name__ == '__main__':
    main()
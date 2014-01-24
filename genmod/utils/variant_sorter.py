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

from genmod.utils.is_number import is_number


class FileSort(object):
    def __init__(self, inFile, outFile=None, sort_mode = 'rank', splitSize=20):
        """ split size (in MB) """
        self._inFile = inFile
        print self._inFile
        
        if outFile is None:
            self._outFile = inFile
        else:
            self._outFile = outFile
                    
        self._splitSize = splitSize * 1000000
        
        self.print_to_screen = True
        
        if sort_mode == 'rank':
            # To sort a CMMS-file on rank score
            self._getKey = lambda variant_line: int(variant_line.rstrip().split('\t')[-1])
        else:
            # to sort a vcf-file on positions
            self._getKey = lambda variant_line: int(variant_line.rstrip().split('\t')[1])
    
    def sort(self):
        files = self._splitFile()

        if files is None:
            """ file size <= self._splitSize """            
            self._sortFile(self._inFile, self._outFile)
            return

        for fn in files:
            self._sortFile(fn)
            
        self._mergeFiles(files)
        self._deleteFiles(files)

        
    def _sortFile(self, fileName, outFile=None):
        lines = open(fileName).readlines()
        get_key = self._getKey
        for line in lines:
            if not is_number(line.rstrip().split('\t')[-1]):
                print line
            elif line.rstrip().split('\t')[-1] == '12':
                print line
            
            
        data = [(get_key(line), line) for line in lines if line!='']
        data.sort(reverse=True)
        lines = [line[1] for line in data]
        if self.print_to_screen:
            print ''.join(lines)
        else:
            if outFile is not None:
                open(outFile, 'w').write(''.join(lines))
            else:
                open(fileName, 'w').write(''.join(lines))
    
    

    def _splitFile(self):
        totalSize = os.path.getsize(self._inFile)
        if totalSize <= self._splitSize:
            # do not split file, the file isn't so big.
            return None

        fileNames = []            
                
        fn,e = os.path.splitext(self._inFile)
        f = open(self._inFile)
        try:
            i = size = 0
            lines = []
            for line in f:
                size += len(line)
                lines.append(line)
                if size >= self._splitSize:
                    i += 1
                    tmpFile = fn + '.%03d' % i
                    fileNames.append(tmpFile)
                    open(tmpFile,'w').write(''.join(lines))
                    del lines[:]
                    size = 0

                                                       
            if size > 0:
                tmpFile = fn + '.%03d' % (i+1)
                fileNames.append(tmpFile)
                open(tmpFile,'w').write(''.join(lines))
                
            return fileNames
        finally:
            f.close()

    def _mergeFiles(self, files):
        files = [open(f) for f in files]
        lines = []
        keys = []
        
        for f in files:
            l = f.readline()        
            lines.append(l)
            keys.append(self._getKey(l))

        buff = []
        buffSize = self._splitSize/2
        append = buff.append
        if not self.print_to_screen:
            output = open(self._outFile,'w')
        try:
            key = max(keys)
            index = keys.index(key)
            get_key = self._getKey
            while 1:
                while key == max(keys):
                    append(lines[index])
                    if len(buff) > buffSize:
                        if self.print_to_screen:
                            print ''.join(buff)
                        else:
                            output.write(''.join(buff))
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
                if self.print_to_screen:
                    print ''.join(buff)
                else:
                    output.write(''.join(buff))
        finally:    
            output.close()

    def _deleteFiles(self, files):   
        for fn in files:
            os.remove(fn)        
    


def main():
    parser = argparse.ArgumentParser(description="Check files.")
    parser.add_argument('infile', type=str, nargs=1, help='Specify the path to the file of interest.')
    parser.add_argument('-out', '--outfile', type=str, nargs=1, help='Specify the path to the outfile.')
    args = parser.parse_args()
    infile = args.infile[0]
    if args.outfile:
        outfile = args.outfile[0]
    else:
        outfile = None
    fs = FileSort(infile, outfile)
    fs.sort()
                    
                
             


if __name__ == '__main__':
    main()
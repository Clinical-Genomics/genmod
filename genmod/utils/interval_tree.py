#!/usr/bin/env python
# encoding: utf-8
"""
interval_tree.py

data: an array of elements where each element contains start coodinate, end coordinate, and element id.
si: index or key of the start coodinate in each element
ei: index or key of the end coordinate in each element
start: position of the start position of the element range
end: posotion of the end position of the element range

for example, a reference genome of a million base pairs with the following features:
features = [[20,400,'id01'],[1020,2400,'id02'],[35891,29949,'id03'],[900000,'id04'],[999000,'id05']]
to make a tree:
myTree = intervalTree(features, 0, 1, 1, 1000000)


Created by unknown.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import pprint

class intervalTree:
    def __init__(self, data, si, ei, start, end):
        '''
        data: an array of elements where each element contains start coodinate, end coordinate, and element id.
        si: index or key of the start coodinate in each element
        ei: index or key of the end coordinate in each element
        start: position of the start position of the element range
        end: posotion of the end position of the element range

        for example, a reference genome of a million base pairs with the following features:
            features = [[20,400,'id01'],[1020,2400,'id02'],[35891,29949,'id03'],[900000,'id04'],[999000,'id05']]
        to make a tree:
            myTree = intervalTree(features, 0, 1, 1, 1000000)
        '''
        self.si = si
        self.ei = ei
        self.start = start
        self.end = end
        self.elementaryIntervals = self.getElementaryIntervals(data, si, ei)
        self.tree = self.recursiveBuildTree(self.elementaryIntervals)
        self.insertData(self.tree, data, si, ei, start, end)
        self.trimTree(self.tree)

    def getElementaryIntervals(self, data, si, ei):
        '''generates a sorted list of elementary intervals'''
        coords = []
        [coords.extend([x[si],x[ei]]) for x in data]
        coords = list(set(coords))
        coords.sort()

        return coords

    def recursiveBuildTree(self, elIntervals):
        '''
        recursively builds a BST based on the elementary intervals.
        each node is an array: [interval value, left descendent nodes, right descendent nodes, [ids]].
        nodes with no descendents have a -1 value in left/right descendent positions.

        for example, a node with two empty descendents:
            [500,                               interval value
                [-1,-1,-1,['id5','id6']],       left descendent
                [-1,-1,-1,['id4']],             right descendent
                ['id1',id2',id3']]              data values

        '''
        center = int(round(len(elIntervals) / 2))

        left = elIntervals[:center]
        right = elIntervals[center + 1:]
        node = elIntervals[center]

        if len(left) > 1:
            left = self.recursiveBuildTree(left)
        elif len(left) == 1:
            left = [left[0],[-1,-1,-1,[]],[-1,-1,-1,[]],[]]
        else:
            left = [-1,-1,-1,[]]

        if len(right) > 1:
            right = self.recursiveBuildTree(right)
        elif len(right) == 1:
            right = [right[0],[-1,-1,-1,[]],[-1,-1,-1,[]],[]]
        else:
            right = [-1,-1,-1,[]]

        return [node, left, right, []]

    def ptWithin(self, pt, subject):
        '''accessory function to check if a point is within a range'''
        if pt >= subject[0] and pt <= subject[1]:
            return True

        return False

    def isWithin(self, query, subject):
        '''accessory function to check if a range is fully within another range'''
        if self.ptWithin(query[0], subject) and self.ptWithin(query[1], subject):
            return True

        return False

    def overlap(self, query, subject):
        '''accessory function to check if two ranges overlap'''
        if self.ptWithin(query[0], subject) or self.ptWithin(query[1], subject) or self.ptWithin(subject[0], query) or self.ptWithin(subject[1], query):
            return True

        return False

    def recursiveInsert(self, node, coord, data, start, end):
        '''recursively inserts id data into nodes'''
        if node[0] != -1:
            left = (start, node[0])
            right = (node[0], end)

            #if left is totally within coord
            if self.isWithin(left, coord):
                node[1][-1].append(data)
            elif self.overlap(left, coord):
                self.recursiveInsert(node[1], coord, data, left[0], left[1])

            if self.isWithin(right, coord):
                node[2][-1].append(data)
            elif self.overlap(right, coord):
                self.recursiveInsert(node[2], coord, data, right[0], right[1])

    def insertData(self, node, data, si, ei, start, end):
        '''loops through all the data and inserts them into the empty tree'''
        for item in data:
            self.recursiveInsert(node, [item[si], item[ei]], item[-1], start, end)

    def trimTree(self, node):
        '''trims the tree for any empty data nodes'''
        dataLen = len(node[-1])

        if node[1] == -1 and node[2] == -1:
            if dataLen == 0:
                return 1
            else:
                return 0
        else:
            if self.trimTree(node[1]) == 1:
                node[1] = -1

            if self.trimTree(node[2]) == 1:
                node[2] = -1

            if node[1] == -1 and node[2] == -1:
                if dataLen == 0:
                    return 1
                else:
                    return 0

    def find(self, node, findRange, start, end):
        '''recursively finds ids within a range'''
        data = []

        left = (start, node[0])
        right = (node[0], end)

        if self.overlap(left, findRange):
            data.extend(node[-1])
            if node[1] != -1:
                data.extend(self.find(node[1], findRange, left[0], left[1]))

        if self.overlap(right, findRange):
            data.extend(node[-1])
            if node[2] != -1:
                data.extend(self.find(node[2], findRange, right[0], right[1]))

        return list(set(data))

    def findRange(self, findRange):
        '''wrapper for find'''
        return self.find(self.tree, findRange, self.start, self.end)

    def pprint(self, ind):
        '''pretty prints the tree with indentation'''
        pp = pprint.PrettyPrinter(indent=ind)
        pp.pprint(self.tree)

def main():
    features = [[20,400,'id01'],[1020,2400,'id02'],[35891,29949,'id03'],[900000,'id04'],[999000,'id05']]
    myTree = intervalTree(features, 0, 1, 1, 1000000)
    print 'Ranges between 200 and 1200:', myTree.findRange([200, 1200])
    print 'Range in only position 34000:', myTree.findRange([900000, 900000])
    print myTree.elementaryIntervals


if __name__ == '__main__':
    main()


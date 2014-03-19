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

class IntervalTree:
    
    def __init__(self, features, start, end):
        """
        data: an array of elements where each element contains start coodinate, end coordinate, and element id.
        start: position of the start position of the element range. Allways 1:based
        end: position of the end position of the element range. Allways 1:based
        
        start index is allways 0 and end index allways 1
        for example, a reference genome of a million base pairs with the following features:
            features = [[20,400,'id01'],[1020,2400,'id02'],[35891,29949,'id03'],[900000,900000,'id04'],[999000,999000,'id05']]
        to make a tree:
            myTree = intervalTree(features, 0, 1, 1, 1000000)
        """
        
        self.start = start
        self.end = end
        self.elementary_intervals = self.get_elementary_intervals(features)
        self.tree = self.recursive_build_tree(self.elementary_intervals)
        self.insert_data(self.tree, features, start, end)
        self.trim_tree(self.tree)

    def get_elementary_intervals(self, features):
        """Generates a sorted list of elementary intervals"""
        coords = []
        try:
            for interval in features:
                if len(interval) != 3:
                    raise SyntaxError('Interval malformed %s. Allways specify start and end position for interval.' % str(interval))
                coords.extend([interval[0],interval[1]])    
        except IndexError:
            raise SyntaxError('Interval malformed %s. Allways specify start and end position for interval.' % str(interval))
        coords = list(set(coords))
        coords.sort()
        return coords

    def recursive_build_tree(self, intervals):
        """
        recursively builds a BST based on the elementary intervals.
        each node is an array: [interval value, left descendent nodes, right descendent nodes, [ids]].
        nodes with no descendents have a -1 value in left/right descendent positions.

        for example, a node with two empty descendents:
            [500,                               interval value
                [-1,-1,-1,['id5','id6']],       left descendent
                [-1,-1,-1,['id4']],             right descendent
                ['id1',id2',id3']]              data values

        """
        center = int(round(len(intervals) / 2))

        left = intervals[:center]
        right = intervals[center + 1:]
        node = intervals[center]

        if len(left) > 1:
            left = self.recursive_build_tree(left)
        elif len(left) == 1:
            left = [left[0],[-1,-1,-1,[]],[-1,-1,-1,[]],[]]
        else:
            left = [-1,-1,-1,[]]

        if len(right) > 1:
            right = self.recursive_build_tree(right)
        elif len(right) == 1:
            right = [right[0],[-1,-1,-1,[]],[-1,-1,-1,[]],[]]
        else:
            right = [-1,-1,-1,[]]

        return [node, left, right, []]

    def pt_within(self, pt, subject):
        """Accessory function to check if a point is within a range"""
        try:
            if pt >= subject[0] and pt <= subject[1]:
                return True
        except TypeError:
            raise TypeError('Interval start and stop has to be integers. %s' % str(subject))

        return False

    def is_within(self, query, subject):
        """Accessory function to check if a range is fully within another range"""
        if self.pt_within(query[0], subject) and self.pt_within(query[1], subject):
            return True

        return False

    def overlap(self, query, subject):
        """Accessory function to check if two ranges overlap"""
        if (self.pt_within(query[0], subject) or self.pt_within(query[1], subject) or 
            self.pt_within(subject[0], query) or self.pt_within(subject[1], query)):
            return True

        return False

    def recursive_insert(self, node, coord, data, start, end):
        """Recursively inserts id data into nodes"""
        if node[0] != -1:
            left = (start, node[0])
            right = (node[0], end)

            #if left is totally within coord
            if self.is_within(left, coord):
                node[1][-1].append(data)
            elif self.overlap(left, coord):
                self.recursive_insert(node[1], coord, data, left[0], left[1])

            if self.is_within(right, coord):
                node[2][-1].append(data)
            elif self.overlap(right, coord):
                self.recursive_insert(node[2], coord, data, right[0], right[1])

    def insert_data(self, node, data, start, end):
        """loops through all the data and inserts them into the empty tree"""
        for item in data:
            self.recursive_insert(node, [item[0], item[1]], item[-1], start, end)

    def trim_tree(self, node):
        """trims the tree for any empty data nodes"""
        data_len = len(node[-1])

        if node[1] == -1 and node[2] == -1:
            if data_len == 0:
                return 1
            else:
                return 0
        else:
            if self.trim_tree(node[1]) == 1:
                node[1] = -1

            if self.trim_tree(node[2]) == 1:
                node[2] = -1

            if node[1] == -1 and node[2] == -1:
                if data_len == 0:
                    return 1
                else:
                    return 0

    def find(self, node, interval, start, end):
        """recursively finds ids within a range"""
        data = []
        
        if len(interval) != 2:
            raise SyntaxError('Interval malformed %s. Allways specify start and end position for interval.' % str(interval))
        
        left = (start, node[0])
        right = (node[0], end)

        if self.overlap(left, interval):
            data.extend(node[-1])
            if node[1] != -1:
                data.extend(self.find(node[1], interval, left[0], left[1]))

        if self.overlap(right, interval):
            data.extend(node[-1])
            if node[2] != -1:
                data.extend(self.find(node[2], interval, right[0], right[1]))

        return list(set(data))

    def find_range(self, interval):
        """wrapper for find"""
        return self.find(self.tree, interval, self.start, self.end)

    def pprint(self, ind):
        """pretty prints the tree with indentation"""
        pp = pprint.PrettyPrinter(indent=ind)
        pp.pprint(self.tree)

def main():
    features = [[20,400,'id01'],[1020,2400,'id02'],[35891,29949,'id03'],[899999,900000,'id04'],[999000,999000,'id05']]
    my_tree = IntervalTree(features, 1, 1000000)
    my_tree.pprint(4)
    print('Ranges between 200 and 1200: %s' % my_tree.find_range([200, 1200]))
    print('Range in only position 90000: %s'  % my_tree.find_range([900000, 900000]))


if __name__ == '__main__':
    main()


#!/usr/bin/env python
'''compare two hhpred outputs'''

import os
import sys
from optparse import OptionParser

def csv_parse(csvfile):
    '''parse the results from the output .hhr file (hhrfile)'''
    file_reader = open(csvfile, "r")
    
    results_dict = {}
    for line in file_reader:
        line = line.strip('\n')
        if line == "":
            continue
        splits = line.split(',')
        
        cluster_id = splits[0]
        hits = int(splits[1])
        results_dict[cluster_id] = hits
       

    file_reader.close()
        
    return results_dict

        
def cross_check(src, dest):
    
    results1 = csv_parse(src)
    results2 = csv_parse(dest)
    both_list = []
    
    for key in results1.keys():
        hits = results1[key]
        if hits > 0:
            continue
        if results2.has_key(key):
            hits2 = results2[key]
            if hits2 > 0:
                continue
        both_list.append(key)
        
    return both_list
         

if __name__ == '__main__':
    usage  = "usage: hhpred_diff.py -s <source file to be compared> -d <dest file to be compared>"
    parser = OptionParser(usage)
    parser.add_option("-s", "--src",  dest="src", type = "string", default=None, help="source file to be compared")
    parser.add_option("-d", "--des",  dest="des", type = "string", default=None, help="dest file to be compared")
    
    (opts, args) = parser.parse_args()
    
    if not (opts.src and os.path.isfile(opts.src) ):
        parser.error("Missing input src file %s"%(opts.src, ))
        sys.exit(1)
        
    if not (opts.des and os.path.isfile(opts.des) ):
        parser.error("Missing input dest file %s"%(opts.des, ))
        sys.exit(1)
        
    both_list = cross_check(opts.src, opts.des)
    
    for item in both_list:
        print "%s,0" % item
        
        
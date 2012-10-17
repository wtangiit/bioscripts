#!/usr/bin/env python
'''compare two hhpred outputs'''

import os
import sys
from optparse import OptionParser

def hhpred_parse(hhrfile):
    '''parse the results from the output .hhr file (hhrfile)'''
    file_reader = open(hhrfile, "r")
    seq_lists = []
    linenumber = 0
    
    results_list = []
    ct = 0
    for line in file_reader:
        line = line.strip('\n')
        if line == "":
            continue
        splits = line.split()
    #    print splits[0]
        if not splits[0].isdigit():
            continue
        ct += 1
        result ={}
        result['Hit'] = splits[1]
        result['Prob'] = float(splits[-9])
        result['E-value'] = splits[-8]
        result['P-value'] = splits[-7]
        result['Score'] = float(splits[-6])
        #print "result=", result
        results_list.append(result)
        if ct == 10:
            break
    file_reader.close()
        
    return results_list
        
def hhpred_diff(src, dest):
    
    results1 = hhpred_parse(src)
    results2 = hhpred_parse(dest)
    
    print "No\tHit\tProb\tE-value\tP-value\tScore"
    
    for i in range(len(results1)):
        print i+1
        print "%s\t%s\t%s\t%s\t%s" % (results1[i]['Hit'],
                                            results1[i]['Prob'],
                                            results1[i]['E-value'],
                                            results1[i]['P-value'],
                                            results1[i]['Score'],)
        print "%s\t%s\t%s\t%s\t%s" % (results2[i]['Hit'],
                                            results2[i]['Prob'],
                                            results2[i]['E-value'],
                                            results2[i]['P-value'],
                                            results2[i]['Score'],)

         

if __name__ == '__main__':
    usage  = "usage: hhpred_diff.py -s <source file to be compared> -d <dest file to be compared>"
    parser = OptionParser(usage)
    parser.add_option("-i", "--inf",  dest="inf", type = "string", default=None, help="input .hhr file")
    parser.add_option("-e", "--evalue",  dest="evalue_thres", type = "float", default=0.001, help="threshold of evalue")
    
    (opts, args) = parser.parse_args()
    
    if not (opts.inf and os.path.isfile(opts.inf) ):
        parser.error("Missing input file %s"%(opts.inf, ))
        sys.exit(1)
        
        
    results_list = hhpred_parse(opts.inf)
    
    hits = 0
    
    for result in results_list:
        evalue = float(result['E-value'])
        if evalue < opts.evalue_thres:
            hits += 1
    print "%s: %s" % (opts.inf, hits)
        
        
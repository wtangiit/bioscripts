#!/usr/bin/env python
'''compare two hhpred outputs'''

import os
import sys
from optparse import OptionParser

def pfam_out_parse(pfam_out_file):
    '''parse the results from pfam_out_file'''
    file_reader = open(pfam_out_file, "r")
    seq_lists = []
    linenumber = 0
    
    results_dict = {}  #{cluster_id:[evalue]}
    
    for line in file_reader:
        line = line.strip()
        if line == "":
            continue
        splits = line.split()
        
        if splits[0] == '#':
            continue
       
        result ={}
        cluster_id = splits[2]
        
        result['Hit'] = splits[0]
        
        try:
            evalue = float(splits[4])
        except ValueError:
            continue
        
        if results_dict.has_key(cluster_id):
            results_dict[cluster_id].append(evalue)
        else:
            results_dict[cluster_id] = [evalue]
        
    file_reader.close()
        
    return results_dict
        
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
    usage  = "usage: hhpred_diff.py -i <input .hhr file> -e <evalue threshold>"
    parser = OptionParser(usage)
    parser.add_option("-i", "--inf",  dest="inf", type = "string", default=None, help="input .hhr file")
    parser.add_option("-e", "--evalue",  dest="evalue_thres", type = "float", default=0.001, help="threshold of evalue")
    
    (opts, args) = parser.parse_args()
    
    if not (opts.inf and os.path.isfile(opts.inf) ):
        parser.error("Missing input file %s"%(opts.inf, ))
        sys.exit(1)
        
    results_dict = pfam_out_parse(opts.inf)
    
    
    for cluster_id in results_dict.keys():
        evalue_list = results_dict[cluster_id]
        hits = 0
        for evalue in evalue_list:
            if evalue < opts.evalue_thres:
                hits += 1
        print "%s,%s" % (cluster_id, hits)
        
        
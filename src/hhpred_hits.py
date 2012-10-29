#!/usr/bin/env python
'''compare two hhpred outputs'''

import os
import sys
from optparse import OptionParser

def find_prob_evalue(str_li):
  #  print str_li
    for i in range(4, len(str_li)):
        try:
            prob = float(str_li[i])
        except ValueError:
            continue
        
        try:
            evalue = float(str_li[i+1])
        except ValueError:
            continue

        return prob, evalue
        

def hhpred_parse(hhrfile):
    '''parse the results from the output .hhr file (hhrfile)'''
    file_reader = open(hhrfile, "r")
    seq_lists = []
    linenumber = 0
    
    results_list = []
    ct = 0
    
    line = file_reader.readline()
    clusterid = line.strip().split()[1]
    for line in file_reader:
        line = line.strip()
        if line == "":
            continue
        splits = line.split()
        
        if not splits[0].isdigit():
            continue
       
        ct += 1
        result ={}
        result['cluster_id'] = clusterid
        result['Hit'] = splits[1]
        
        prob, evalue = find_prob_evalue(splits)
        result['Prob'] = prob
        result['E-value'] = evalue
        
      #  print "hit=%s, prob=%s, e-value=%s" % (result['Hit'], prob, evalue)
            
        results_list.append(result)
        #print result
        if ct == 20:
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
    usage  = "usage: hhpred_diff.py -i <input .hhr file> -e <evalue threshold>"
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
       # print "evalue=", evalue
        if evalue < opts.evalue_thres:
            hits += 1
    print "%s,%s" % (result['cluster_id'], hits)
        
        
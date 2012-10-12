#!/usr/bin/env python

import sys, os, random, re
from optparse import OptionParser
from Bio import SeqIO

if __name__ == '__main__':
    usage  = "usage: %prog -i <file1> -o <rejectfile> \n paired-to-interleaved, outputs to std out"
    parser = OptionParser(usage)
    parser.add_option("-i", "--in",  dest="inp", default=None, help="input file")
    parser.add_option("-o", "--out",  dest="outp", default=None, help="output file")
    parser.add_option("-n", "--num",  dest="num", default=1 , help="number of columns to keep in header")
    
    (opts, args) = parser.parse_args()
    
    if opts.num:
        num = int(opts.num)
    else:
        num = 1
        
    if not (opts.inp and os.path.isfile(opts.inp) ):
        parser.error("Missing input file" )
        
    in_fh  = open(opts.inp)
    
    if(opts.outp):
        out_fh= open(opts.outp,"w")
    else:
        out_fh = open(opts.inp+".new", "w") 
    
    for line in in_fh:
        if line[0] == '>':
            segs = line.split(' ')
            new_segs = segs[0:num]
            new_header= " ".join(new_segs)
            new_header+='\n'
            out_fh.write(new_header)
        else:
            out_fh.write(line)
            
    in_fh.close()
    out_fh.close()

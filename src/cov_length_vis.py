#!/usr/bin/env python

import sys, os, random, re
from optparse import OptionParser
from Bio import SeqIO
import matplotlib
import matplotlib.pyplot as plt

if __name__ == '__main__':
    usage  = "usage: %prog -i <inputfile>"
    parser = OptionParser(usage)
    parser.add_option("-i", "--in",  dest="inp", default=None, help="input file")
#    parser.add_option("-o", "--out",  dest="outp", default=None, help="output file")
  
    (opts, args) = parser.parse_args()
        
    if not (opts.inp and os.path.isfile(opts.inp) ):
        parser.error("Missing input file" )
        
    in_fh  = open(opts.inp)
        
    specs = []
    
    records=SeqIO.parse(in_fh, "fasta")
    for seq_record in records:
        spec = {}
        slices = seq_record.description.split(" ")
        spec['id'] = slices[0]
        for i in range(1, len(slices)):
            strs = slices[i].split("=")
            spec[strs[0]] = strs[1]
        specs.append(spec)
        
    def cmp_cov(spec1, spec2):
        return -cmp(float(spec1["avgcov"]), float(spec2["avgcov"]))
    
    specs.sort(cmp_cov)
    
    length_list = []
    
    acc_lengths = []
    acc_length = 0
    for spec in specs:
        length = int(spec["length"])
        acc_length += length
        acc_lengths.append(acc_length) 
    
    x = range(len(acc_lengths))
    y = acc_lengths
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.set_xlabel("coverage rank")
    ax.set_ylabel("accumulate length")
    ax.set_title(opts.inp)
    plt.savefig(opts.inp+".png")
        
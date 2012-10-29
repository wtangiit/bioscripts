#!/usr/bin/env python

import sys, os, random, re
from optparse import OptionParser
from Bio import SeqIO
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

GRID_LINEWIDTH = 0.6
GRID_ALPHA = 0.3
PLOTWIDTH = 1.4
    
def plot_scatter(specs, x_key, y_key, name):
        
    x_list = []
    y_list = []
    size_list = []    

    for spec in specs:
        
        x_value = int(float(spec[x_key]))
        y_value = int(float(spec[y_key]))
        
        x_list.append(x_value)
        y_list.append(y_value)
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
#    ax.set_yscale('log')
    ax.plot(x_list, y_list, '.', linewidth=0.5)
    ax.set_xlabel(x_key)
    ax.set_ylabel(y_key)
    ax.set_title(name)
#    plt.show()
    plt.savefig("scatter_%s-%s-%s.png" % (name, y_key, x_key))
    
    
def plot_hist(specs, key, name):
        
    x_list = []
    for spec in specs:
        x_value = int(float(spec[key]))
        #x_list.append(x_value)
        x_list.append(math.log(x_value, 10))
        #print x_value
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
#    ax.set_yscale('log')

    

    ax.hist(x_list, bins=10, facecolor='green')
    
    ax.set_xlabel("log(%s)" % key)
    ax.set_title(name)
   # plt.show()
    plt.savefig("hist_%s-log(%s).png" % (name, key))
    

  # plt.savefig("scatter_%s-%s-%s-log.png" % (name, y_key, x_key))

def plot_acc_length(specs, name):
        
    def cmp_cov(spec1, spec2):
        if spec1.has_key("cov"):
            return -cmp(float(spec1["cov"]), float(spec2["cov"]))
        else:
            return 0
        
    def cmp_length(spec1, spec2):
        return -cmp(float(spec1["length"]), float(spec2["length"]))
    
    specs.sort(cmp_cov)
    #specs.sort(cmp_length)
    
    length_list = []
    cov_list = []    
    
    acc_lengths = []
    length = 0
    for spec in specs:
        length = int(spec["length"])
        acc_length += length
        acc_lengths.append(acc_length)
        
        if spec.has_key("cov"):
            cov = spec["cov"]

            
        cov_list.append(cov)      
        length_list.append(length)  
    
    x = range(len(acc_lengths))
    y = acc_lengths
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.set_xlabel("coverage rank")
    ax.set_ylabel("accumulate length")
    ax.set_title(name)
    plt.savefig(name+".png")

    
    
def plot_acc_lencov(specs, name):
        
    def cmp_cov(spec1, spec2):
        if spec1.has_key("cov"):
            return -cmp(float(spec1["cov"]), float(spec2["cov"]))
        else:
            return 0
        
    def cmp_length(spec1, spec2):
        return -cmp(float(spec1["length"]), float(spec2["length"]))
    
    specs.sort(cmp_cov)
    #specs.sort(cmp_length)
    
    length_list = []
    cov_list = []    
    
    acc_lengths = []
    acc_length = 0
   
    lencov_list = []
    
    acc_lencov = 0
    
    #calc total len*cov
    for spec in specs:
        length = int(spec["length"])
        acc_length += length
        acc_lengths.append(acc_length)

        cov = float(spec["cov"])
            
        cov_list.append(cov)      
        length_list.append(length)
        
        lencov = length*cov
        size = int(spec["size"])
        #print lencov / size 
        lencov_list.append(lencov)
                
    acc_fexplain_list = []
    acc_fexplain = 0
    sum_lencov = sum(lencov_list)
    for item in lencov_list:
        fexplain = item / sum_lencov
        acc_fexplain += fexplain
        acc_fexplain_list.append(acc_fexplain)
            
    x = range(len(cov_list))
    y = acc_fexplain_list
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.set_xlabel("coverage rank")
    ax.set_ylabel("accumulate fexplain [ len*cov / sum(len*cov) ]")
    ax.set_title(name)
    ax.yaxis.grid(True, linestyle='-', which='major', alpha=GRID_ALPHA)
#    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.grid(True, linestyle='-', which='major', alpha=GRID_ALPHA)
    plt.savefig(name+"acc_fexplain.png")

if __name__ == '__main__':
    usage  = "usage: %prog -i <inputfile>"
    parser = OptionParser(usage)
    parser.add_option("-i", "--in",  dest="inp", default=None, help="input file")
    parser.add_option("-c", "--cdf", dest="cdf", default=None, help="draw cdf")
    parser.add_option("-x", "--x",  dest="x", default=None, help="name for x_axis")
    parser.add_option("-y", "--y",  dest="y", default=None, help="name for y_axis")
    parser.add_option("-t","--hist",  dest="hist", default=None, help="histogram")
    
    parser.add_option("-a", "--acc", dest="acc", action="store_true", default=False, help="draw accumulating function")   
  
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
            key = strs[0]
            value = strs[1]
            if key=="avgcov":
                key = "cov"
            spec[key] = value
        specs.append(spec)
    if (opts.x and opts.y):
        plot_scatter(specs, opts.x, opts.y, opts.inp)
        #plot_hist(specs, opts.x, opts.y, opts.inp)
        
    if (opts.hist):
        plot_hist(specs, opts.hist, opts.inp)
        
    if opts.acc:
        plot_acc_lencov(specs, opts.inp)
        #plot_acc_length(specs, opts.inp)
           
#!/usr/bin/env python

'''
calc.pwm.py
calculate PWM scores for given sequences and start/stop matrix
''' 


import os
import math
from optparse import OptionParser

nt_list = ['A', 'C', 'G', 'T']

nt_dict = {'A':0, 'C': 1, 'G':2, 'T':3, 
           'a':0, 'c':1, 'g':2, 't':3}

complement_dict = {'A':'T', 'C': 'G', 'G':'C', 'T':'A', 
           'a':'t', 'c':'g', 'g':'c', 't':'a'}

STRATIFY = False

MIN_GC_CONTENT = 26
MAX_GC_CONTENT = 70
NUM_STRATIFY = 45
NUM_M_STATE = 6
NUM_DIMER = 16
NUM_NT = 4

def trimer_to_int(triplet):
    '''return number by triplet'''
    t1 = nt_dict.get(triplet[0])
    t2 = nt_dict.get(triplet[1])
    t3 = nt_dict.get(triplet[2])
    
    if t1 >= 0 and t2 >=0 and t3 >=0:
        return t1 * 16 + t2 * 4 + t3
    else:
        return -1
    
def get_gc_content(sequence):
    '''return gc_content% of a given sequence'''
    gc_count = 0
    for ch in sequence:
        if ch in ['G', 'C', 'g', 'c']:
            gc_count += 1
#    print "%f\t%f\t%f"%(float(gc_count) / len(sequence), round(float(gc_count) / len(sequence), 2) * 100, int(round(float(gc_count) / len(sequence), 2) * 100)  ) 
    gc_content = int(float(gc_count) / len(sequence) * 100 + 0.5)
    #print "gc countent %s rounded to %s" % (float(gc_count) / len(sequence) * 100, gc_content)
    if gc_content < MIN_GC_CONTENT:
        gc_content = MIN_GC_CONTENT
    if gc_content > MAX_GC_CONTENT:
        gc_content = MAX_GC_CONTENT
    return gc_content

def parse_input_file(filename, noncoding=False):
    '''parse the input file'''
    infile = open(filename, "r")
    seq_lists = []
    linenumber = 0    
    for line in infile:
        line = line.strip('\n')
        line = line.strip('\r')
        splits = line.split('\t')
        linenumber += 1
        if len(splits) == 3:
            seq_len = len(splits[2])
            if seq_len > 0:
                if not noncoding and seq_len < 123:
                    print "**********************************************"
                    print "Warning!!:the input data contains invalid data, line %d, sequence length=%s" % (linenumber, seq_len)
                    print "The invalid sequence is thrown out to continue, but replacing input data and re-training is suggested." 
                    print "**********************************************"
                    continue
                seq_lists.append(splits[2])
    return seq_lists
            
def get_reverse_complement(seq):
    '''return the reverse complement of the given sequence'''
    seq = seq[::-1]
    rseq= ""
    for ch in seq:
        try:
            rseq += complement_dict[ch]
        except KeyError:
            rseq += "N"
    return rseq            

def parse_matrix(matrix_file):
    '''parse matrix file  '''
    filehandle = open(matrix_file, "r")
    pos_index = 0  # 0 to 60
    
    matrix = [[[0 for t in range(64)] for p in range(61)] for g in range(45)]
                         
    for line in filehandle:
        line = line.strip('\n')
        line = line.strip('\r')
        numbers = line.split()
        if len(numbers) == 1:
            gc_index = int(numbers[0]) - 26
            pos_index = 0
        else:
            for triplet_index in range(64):
                matrix[gc_index][pos_index][triplet_index] = numbers[triplet_index]
            pos_index += 1
    return matrix

def calc_pwm_score(seq_list, matrix, type="start"):
    '''calculate pwm scores, return a list of all the scores'''
        
    score_dict = {}    
    
    for seq in seq_list:
        
        gc_content = get_gc_content(seq)
        
        if not score_dict.has_key(gc_content):
            score_dict[gc_content] = []
        
        if type=="start":
            subseq = seq[30:93]
        elif type=="end":
            subseq = seq[-123:-60]
        elif type=="start1":
            rc_seq = get_reverse_complement(seq)
            subseq = rc_seq[-93:-30]
        elif type=="end1":
            rc_seq = get_reverse_complement(seq)
            subseq = rc_seq[60:123]
        else:
            subseq = seq
        
        score = 0
        for i in range(61):
            
            #skip stop codon when calculating (start codon is not skipped)
            if type=="end" and i==60:
                continue
            elif type=="start1" and i==0:
                continue
                
            trimer = subseq[i:i+3]
            trimer_index = trimer_to_int(trimer)
            score -= math.log(float(matrix[gc_content-26][i][trimer_index]))
        
        score_dict[gc_content].append(score)
    return score_dict
                                    
            
if __name__ == '__main__':
    usage  = "usage: %prog -i <input sequence file> -n <input noncoding file> [-g]"
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",  dest="input", type = "string", default=None, help="Input gene sequence file.")
    parser.add_option("-m", "--matrix",  dest="matrix", type = "string", default=None, help="start or stop matrix.")
    parser.add_option("-t", "--type",  dest="type", type = "string", default=None, help="type of matrix (start, stop, start1, stop1).")
    parser.add_option("-g", "--gc", dest="gc_content", action="store_true", default=False, help="stratify by gene GC content")
    
    (opts, args) = parser.parse_args()
    
    STRATIFY = opts.gc_content
    
    msg = "Stratify= %s" % STRATIFY
    if STRATIFY==False:
        msg += ", use -g to enable gc_content stratification"
    print msg    
 
    if not (opts.input and os.path.isfile(opts.input) ):
        parser.error("Missing input file %s"%(opts.input, ))
        
    inputFile = opts.input 
    
    seq_list = parse_input_file(inputFile)
    print "total # of sequences=", len(seq_list)
    
    
    if opts.matrix:
        matrix = parse_matrix(opts.matrix)
        score_dict = calc_pwm_score(seq_list, matrix)
        
        for i in range(26, 71):
            print i
            if score_dict.has_key(i):
                print score_dict[i]

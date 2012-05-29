#!/usr/bin/env python

'''inject substitution error to amino acid sequences at given error rate (%)'''

import sys
import random

__helpmsg__="Usage: error_inject.py input_file.fa error_rate"

error_types = ["sub", "insert", "delete"]

aa_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
'P', 'S', 'T', 'W', 'Y', 'V']

def aa_substitute(old_aa):
    '''give an old_aa, subtitute a new aa randomly'''
    
    ran = int(random.random() * 100) % 20
    while old_aa == aa_list[ran]:
        ran = int(random.random() * 100) % 20
    return aa_list[ran]
        
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print __helpmsg__
        sys.exit(1)
    args = sys.argv
    
    input_file = args[1]
    error_rate = float(args[2])
        
    print "error_rate = ", error_rate
        
    infile = open(input_file, "r")
    outfile_name = "%s.%s" % (input_file, args[2])
    outfile = open(outfile_name, "w")
    
    for line in infile:
        line = line.strip('\n')
        line = line.strip('\r')
        if line[0] == '>':
            outfile.write(line + '\n')
            continue
        
        newline = ""
        for residu in line:
            ran = random.random()
            if ran < error_rate:
                new_posistion = aa_substitute(residu)
            else:
                newline += residu
        #print newline
        outfile.write(newline + '\n')
    outfile.close()        

#!/usr/bin/env python

import sys, os, string
from optparse import OptionParser
from Bio import SeqIO


if __name__ == '__main__':
  usage  = "usage: %prog -i <input sequence file> -t <table file>\nfunction: splits similairities by genome in fasta header"
  parser = OptionParser(usage)
  parser.add_option("-i", "--input",  dest="input", default=None, help="Input sequence file.")
  parser.add_option("-t", "--table",  dest="table", default=None, help="table of identifiers")
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
  
  (opts, args) = parser.parse_args()
  if not (opts.input and os.path.isfile(opts.input) ):
    parser.error("Missing input file" )

  if opts.verbose: sys.stdout.write("Processin %s ... "%opts.input)
  in_handle  = open(opts.input)
  
  in_handle   = open(opts.input, "r")
  tbl_handle  = open(opts.table, "r")
  out_handles={}
  for line in tbl_handle:
     a=opts.input +"."+ line.rstrip()
     print "Opening ",a
     out_handles[line.rstrip()] = open(a , "w")

  something=SeqIO.parse(in_handle, "fasta") 
  for seqrecord in something: 
    try:
#    try out_handles[seqrecord.description[0:9]] > 0 :
      indx = string.find(str(seqrecord.description), "NC_") 
      if indx >=  0 :
        acc = str(seqrecord.description[indx:indx+9])
#        print "acc", acc , "idx ", indx
      out_handles[acc].write(">%s\n%s\n"%(seqrecord.description, str(seqrecord.seq)) )
#      sys.stdout.write(str(seqrecord.description)[0:9]+"\n")
#      sys.stdout.write(">%s\n%s\n"%(seqrecord.description, str(seqrecord.seq)) )
    except: 
        pass
	

#!/usr/bin/env python

import sys, os
from optparse import OptionParser
from Bio import SeqIO

if __name__ == '__main__':
  usage  = '''usage: %prog -i <input sequence file> -p <input ptt> [-f] [-r] \n    p
urpose: produces tables of sequence subsets corresponding to genes '''
  parser = OptionParser(usage)
  parser.add_option("-i", "--input",   dest="input", default=None, help="Input sequence file.")
  parser.add_option("-p", "--ptt",     dest="ptt", default=None, help="Input ptt table.")
  parser.add_option("-b", "--buffer",  dest="buf", default=60, help="Forward / reverse buffer region")
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
  parser.add_option("-f", "--fasta", dest="fasta", action="store_true", default=False, help="Fasta output (default csv)")
  parser.add_option("-t", "--tab", dest="tabsep", action="store_true", default=False,  help="Tab separated")

  (opts, args) = parser.parse_args()
  buf = int(opts.buf)
  if not (opts.input and os.path.isfile(opts.input) ):
    parser.error("Missing input file %s"%(opts.input, ))
  if not (opts.ptt and os.path.isfile(opts.ptt) ):
    parser.error("Missing input file %s"%(opts.ptt, ))
  upstream = buf
  downstream = buf
  if opts.verbose: sys.stderr.write("Processing %s and %s... \n"%(opts.input, opts.ptt))
  in_handle  = open(opts.input)
  ptt_handle = open(opts.ptt)
  record=SeqIO.parse(in_handle, "fasta").next()
  for line in ptt_handle:
     fields = line.split()
     field1 = fields[0]
     field2 = field1.split(".")
     if len(field2) == 3 : 
       start = int(field2[0])
       stop  = int(field2[2])
       direction = fields[1]
       try:
         label= fields[5]
       except:
         label=""     
       if direction == "+":
        if start-1-upstream < 0 or stop+downstream > len(record.seq) :
          sys.stderr.write("problem with %s out of range. \n"%(label))    
        else:
          if opts.fasta:
            print ">%s_%s-%s_%s_plusminus%d length=%d"%(label, start, stop, direction, buf, (stop-start))
          else:
            print "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start)),
          seq = str(record.seq[(int(start)-1-upstream):(int(stop)+downstream)])
          print seq 
       if direction == "-":
        if start-1-downstream  < 0 or stop + upstream > len(record.seq)  :
          sys.stderr.write("problem with %s out of range. \n"%(label))    
        else: 
          if opts.fasta:
            print ">%s_%s-%s_%s_plusminus%d length=%d"%(label, start, stop, direction, buf, (stop-start))
          else:
            print "%s_%s-%s_%s_plusminus%d\t%d\t"%(label, start, stop, direction, buf, (stop-start)),
          seq = str(record.seq[int(start)-1-downstream:int(stop)+upstream].reverse_complement() )
          print seq
  in_handle.close()

  if opts.verbose: sys.stderr.write("Done. \n")

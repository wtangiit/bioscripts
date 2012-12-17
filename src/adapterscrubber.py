#!/usr/bin/env python

# performs exhaustive ungapped alignment against a database of "adapter" sequences

import sys, os,string
from optparse import OptionParser
from Bio import SeqIO

removeambigtable= string.maketrans('RYWSMKHBVDNrywsmkhbvdn','N'* 22)

def ip (a,b, d):   # takes two sequences and offset and returns (number of matches, overlap length)
   assert len(b) >= len(a)
   if d < 0:
      r = len(a) + d 
   elif (len(a) + d) > len(b):
      r =  len(b) - d
   else: 
      r= len(a) 
   c=0
   for i in range( 0, r ) :
      if (a[i] == b[i+d] and a[i] != "N") :
         c+=1
   if r==0:  r=1
   return c ,r

def align ( sequence, adapter ):  # takes two sequences, returns best alignment (number of matches, alignlength, offset)
   if len(sequence) > len(adapter) :
      a, b = adapter, sequence   
   else:
      a, b = sequence, adapter
   assert len(b) >= len(a)
   la = len(a)
   lb = len(b)
   bestm = 0
   bestr = 0
   besto = 0
   for i in range(-len(a), len(b)):
      (m,r) = ip(a, b, i ) 
      c = float(m)/float(r)
      if m > bestm and c > MINALIGNID and r > MINOVERLAP :
         bestm = m
         besto = i
         bestr = r
   return (bestm , bestr, besto)

def bestalign( sequence, adapters):
   type(sequence)
   type(adapters)
   besta = 0
   bestk = 0
   for key in adapters.keys():
      (m,r,o) = align(sequence, adapters[key]) 
      if m > besta:
         besta = m
         bestk = key 
   return (besta, bestk)

if __name__ == '__main__':
  usage  = "usage: %prog -d <database sequence file> <input fatsta file>"
  parser = OptionParser(usage)
  parser.add_option("-i", "--input",  dest="input", default=None, help="Input sequence file.")
  parser.add_option("-d", "--database", dest="database", default=None, help="Database fasta of adapter sequences")
  parser.add_option("-m", "--minoverlap", dest="MINOVERLAP", default=10, help="Minimum overlap paramter")
  parser.add_option("-f", "--fractionid", dest="MINALIGNID", default=.9, help="Minimum alignment id")
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
  
  (opts, args) = parser.parse_args()
  try:
    filename = args[0]
  except:
    parser.error("Input file must be specified")

  MINOVERLAP=int(opts.MINOVERLAP)
  MINALIGNID=float(opts.MINALIGNID)
  dbfile = opts.database

  if not (filename and os.path.isfile(filename) ):
    parser.error("Missing input file" )
  if not (dbfile and os.path.isfile(dbfile) ):
    parser.error("Missing input file" )

  if opts.verbose: sys.stdout.write("Testing %s against adapter database %s\n"%(filename, dbfile))
  adapters = {}
  for seq_record in SeqIO.parse(dbfile , "fasta"):
     adapters[seq_record.description] =        string.upper(str(seq_record.seq).translate(removeambigtable))
     adapters["%s.R"%seq_record.description] = string.upper(str(seq_record.seq.reverse_complement()).translate(removeambigtable))

  for seq_record in SeqIO.parse(filename , "fasta"):
     (a,b) = bestalign( seq_record.seq, adapters )

     flag = a > 0
#     print seq_record.description, a, len(adapters[b]), b
     print seq_record.description, str(int(flag))




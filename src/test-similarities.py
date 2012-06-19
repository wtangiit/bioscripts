#!/usr/bin/env python

import sys, os
from optparse import OptionParser
from Bio import SeqIO
dbpath="../"

if __name__ == '__main__':
  usage  = "usage: %prog -s <stem> -t <table>\nfunction: runs similarities against correct answer "
  parser = OptionParser(usage)
#  parser.add_option("-i", "--input",  dest="input", default=None, help="Input sequence file.")
  parser.add_option("-s", "--stem",  dest="stem", default=None, help="input file Stem")
  parser.add_option("-t", "--table",  dest="table", default=None, help="table of identifiers")
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
  
  (opts, args) = parser.parse_args()
  if not (opts.table and os.path.isfile(opts.table) ):
    parser.error("Missing table file" )

  if opts.verbose: sys.stdout.write("Processing %s ... \n"%opts.table)
  tbl_handle  = open(opts.table, "r")
  out_handles={}
  stem=opts.stem

  for line in tbl_handle:
     accession = line.rstrip()
     a=opts.stem +"."+ line.rstrip()
     if os.path.isfile(dbpath+accession+".faa.phr") != True:
       print "Can't find "+ dbpath+accession+".faa.phr"
       os.system("formatdb -i %s%s.faa"%(dbpath, accession)) 
     else:
       print "Found "+ dbpath+accession+".faa.phr"

#     print "Opening ",a
     if os.path.isfile("%s.%s.out"%(stem, accession)) != True:
       print("blastall -p blastp -F F -d %s%s.faa -i %s.%s -o %s.%s.out -m 8 -a 8 "% (dbpath,accession, stem, accession, stem, accession)) 
       os.system("blastall -p blastp -F F -d %s%s.faa -i %s.%s -o %s.%s.out -m 8 -a 8 "% (dbpath,accession, stem, accession, stem, accession)) 
     else:
       print "File %s.%s.out already exists"%(stem, accession)



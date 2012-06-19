#!/usr/bin/env python

import sys, os
from optparse import OptionParser
from Bio import SeqIO
dbpath="../"

if __name__ == '__main__':
  usage  = "usage: %prog -s <stem> -t <table>\nfunction: adds up similarities tables for accuracy"
  parser = OptionParser(usage)
#  parser.add_option("-i", "--input",  dest="input", default=None, help="Input sequence file.")
  parser.add_option("-s", "--stem",  dest="stem", default=None, help="input file Stem")
  parser.add_option("-t", "--table",  dest="table", default=None, help="table of identifiers")
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
  
  (opts, args) = parser.parse_args()
  if not (opts.table and os.path.isfile(opts.table) ):
    parser.error("Missing table file" )

  if opts.verbose: sys.stdout.write("Processing %s ... "%opts.table)
  tbl_handle  = open(opts.table, "r")
  out_handles={}
  stem=opts.stem
  print "cat %s | fastalengths.pl | sumit.pl "%(stem,)
  countaa = os.popen("cat %s | fastalengths.pl | sumit.pl "%(stem,)).read().rstrip()
  countsq = os.popen("grep -c '>' %s "%(stem,)).read().rstrip()
#  print "Countaa : ", countaa, "countsq: ", countsq
  print "accession aa-pred sq-pred aa-anno sq-anno lablel  aa-uniq-pred sq-uniq-pred aa-uniq-anno sq-uniq-pred"
  for line in tbl_handle:
     accession = line.rstrip()
     a=opts.stem +"."+ accession 
     out=opts.stem +"."+ accession +".out"
     correctfile = dbpath+accession + ".faa"
#     print "Opening ",a
     countaaL = os.popen("cat %s | fastalengths.pl | sumit.pl "%(a,)).read().rstrip()
     countsqL = os.popen("grep -c '>' %s "%(a,)).read().rstrip()
     countaaC = os.popen("cat %s | fastalengths.pl | sumit.pl "%(correctfile,)).read().rstrip()
     countsqC = os.popen("grep -c '>' %s "%(correctfile,)).read().rstrip()
     uniq1  = os.popen("cat %s | sort  -k 1,1 -k 4,4rn | uniq-col1.pl | cut -f 4 | sumit.pl"%(out,) ).read().rstrip()
     uniq1C = os.popen("cat %s | sort  -k 1,1 -k 4,4rn | uniq-col1.pl | wc -l "%(out,) ).read().rstrip()
     uniq2  = os.popen("cat %s | sort  -k 2,2 -k 4,4rn | uniq-col2.pl | cut -f 4 | sumit.pl"%(out,) ).read().rstrip() 
     uniq2C = os.popen("cat %s | sort  -k 2,2 -k 4,4rn | uniq-col2.pl | wc -l "%(out,) ).read().rstrip() 
     aaSn  = float(uniq2) / float(countaaC)
     aaPPV = float(uniq1) / float(countaaL)
#     print accession, "Countaa : ", countaaL, "countsq: ", countsqL, "out ", out, "un1 ", uniq1, "un2 ", uniq2 
     print accession, countaaL, countsqL, countaaC, countsqC, out,  uniq1, uniq1C, uniq2 , uniq2C, "%.04f"%(aaSn,), "%.04f"%(aaPPV,)
#    print "accession aa-pred sq-pred aa-anno sq-anno lablel  aa-uniq-pred sq-uniq-pred aa-uniq-anno sq-uniq-pred"

#     print("blastall -p blastp -F F -d %s%s.faa -i %s.%s -o %s.%s.out -m 8 -a 8 "% (dbpath,accession, stem, accession, stem, accession)) 
#     print("blastall -p blastp -F F -d %s%s.faa -i %s.%s -o %s.%s.out -m 8 -a 8 "% (dbpath,accession, stem, accession, stem, accession)) 

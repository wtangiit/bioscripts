#!/usr/bin/env python
# General-purpose get-the-data script; tries to recognize accesison number and download data accordingly.  

# For GENBANK, depends on Biopython  example:   fetch.py NC_000907 -f gb
# For MGRAST and WGS depends on curl example:   fetch.py AAAE 
#                                               fetch.py 4440055.3 
# For SRA depends on aspera ascp in ASPERAPATH  fetch.py SRR000311 

import sys, os, re
from optparse import OptionParser


def retrieveGENBANKbyaccession(accession, format="fasta"):
  from Bio import Entrez
  sys.stderr.write("Downloading %s, requested format %s\n"%(accession, format))
  Entrez.email = "trimble@anl.gov"
  handle = Entrez.efetch(db="nucleotide", id=accession, rettype=format)
  if format=="fasta":
    f = open("%s.fna"% accession, "w")
  elif format=="gb":
    f = open("%s.gbk"% accession, "w")
  else :
    sys.stderr.write("Warning: unrecognized format %s, defaulting to XML\n"%format)
    f = open("%s.xml"% accession, "w")
  f.write(handle.read())
  f.close()

def retrieveWGSbyaccession(accession, format="fasta"):
  print "four-digit accession %s"%accession
  s="curl ftp://ftp.ncbi.nih.gov/genbank/wgs/wgs.%s.1.fsa_nt.gz >  %s.fna.gz"%(accession, accession)
  os.system(s) 

def retrieveMGRbyaccession(accession, format="fasta"):
#  http://api.metagenomics.anl.gov/reads/mgm4447971.3
  a=re.search("^(4......\..)$", accession).group(1)
  if key=="":
    sys.stderr.write("MGR webkey not defined\n")
    s="curl http://api.metagenomics.anl.gov/reads/mgm%s/ -D /tmp/fetch-dump > %s.gz"%(a, a) 
  else: 
    sys.stderr.write("Using MGR webkey %s\n"%key)
    s="curl 'http://api.metagenomics.anl.gov/reads/mgm%s/?&auth=%s' -D /tmp/fetch-dump > %s.gz"%(a, key, a)
  sys.stderr.write("Executing %s\n"%s) 
  print s
  os.system(s )
  
def retrieveSRRbyaccession(accession, format="fastq"):
  try:
    asperapath=os.environ["ASPERAPATH"]
  except KeyError:
    sys.exit("Error: environment variable ASPERAPATH must be defined and point to aspera to use SRR download")
  a= re.search(".RR(......)", accession).group(1)


  print  "Fetching SRR accession %s"%accession
  tla      = accession[0:3]
  stem     = accession[0:6]
  filename = accession+".sra"
  s = "%s/bin/ascp -i %s/etc/asperaweb_id_dsa.putty anonftp\@ftp-trace.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s ./%s" % ( asperapath, asperapath, tla, stem, accession, filename, filename)
  print "with %s"%(s)
  if int(os.popen("ascp 2>&1 |wc").read().split()[0]) > 10:
    os.system(s)
  else: 
    sys.exit("Sorry, can't find ascp.\t")
  print s 
  s = "fastq-dump --split-3 ./%s "%(filename)
  os.system(s)

if __name__ == '__main__':
  usage  = "usage: %prog <accession number> [-f <format>]"
  parser = OptionParser(usage)
  parser.add_option("-f", "--format", dest="format", default=None, help="Data format")
  (opts, args) = parser.parse_args()
  try:
    key=os.environ["MGRKEY"]
  except KeyError:
    key=""

  try :
    accession = args[0]
  except IndexError:
    parser.error.exit("accession is a required parameter\n%s"%usage)

  if re.search("^(4......\..)$", accession):  # MGR accession
    retrieveMGRbyaccession(accession)
  elif re.search(".RR(......)", accession):   # SRR accession
    retrieveSRRbyaccession(accession)
  elif re.search("^(\w\w\w\w)$", accession):  # WGS accession
    retrieveWGSbyaccession(accession)
  elif re.search("^(N._\d\d\d\d\d\d)$", accession):  # WGS accession
    retrieveGENBANKbyaccession(accession, format=opts.format)
  else: 
    print "Don't recognize acecssion %s"%accession
    sys.exit()

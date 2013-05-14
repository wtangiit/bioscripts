#!/usr/bin/env python
# General-purpose get-the-data script; tries to recognize accesison number and download data accordingly.  

# For GENBANK, depends on Biopython  example:   fetch.py NC_000907 -f gb
# For MGRAST and WGS depends on curl example:   fetch.py AAAE 
#                                               fetch.py 4440055.3 
# For SRA depends on aspera ascp in ASPERAPATH  fetch.py SRR000311 

#ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP009/SRP009660 

#ftp://ftp.sra.ebi.ac.uk/vol1/ERA007/ERA007448/srf   # example
#ftp://ftp.sra.ebi.ac.uk/vol1/SRR387/SRR387449   #  fails
#ftp://ftp.sra.ebi.ac.uk/vol1/SRA048/SRA048196   #fails

#ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR387/SRR387449  # seems to return gz
# per instructions at http://www.ebi.ac.uk/ena/about/sra_data_download
# ascp -QT -i /homes/trimble/build/aspera/etc/asperaweb_id_dsa.putty era-fasp@fasp.sra.ebi.ac.uk:vol1/srr/SRR387/SRR387449 ./SRR387449-2.gz
# works as a fallback!

import sys, os, re
from optparse import OptionParser


def retrieveGENBANKbyaccession(accession, rformat="fasta"):
    from Bio import Entrez
    if rformat == None:
        rformat = "fasta"
    sys.stderr.write("Downloading %s, requested format %s\n" % (accession, rformat))
    Entrez.email = "trimble@anl.gov"
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype=rformat)
    if rformat == "fasta":
        f = open("%s.fna" % accession, "w")
    elif rformat == "gb" or rformat == "gbk" or rformat == "genbank":
        f = open("%s.gbk" % accession, "w")
    elif rformat == "ft":
        f = open("%s.ft" % accession, "w")
    else :
        sys.stderr.write("Warning: unrecognized format %s, defaulting to XML\n" % rformat)
        f = open("%s.xml"% accession, "w")
    f.write(handle.read())
    f.close()

def retrieveWGSbyaccession(accession, rformat="fasta"):
    print "four-digit accession %s" % accession
    s = "curl ftp://ftp.ncbi.nih.gov/genbank/wgs/wgs.%s.1.fsa_nt.gz >  %s.fna.gz" % (accession, accession)
    os.system(s) 

def retrieveMGRbyaccession(accession, rformat="fasta"):
#  http://api.metagenomics.anl.gov/sequences/mgm4447971.3 OBSOLETE
#  http://api.metagenomics.anl.gov/reads/mgm4447971.3   OBSOLETE
    a = re.search("^(4......\..)$", accession).group(1)
    if key == "":
        sys.stderr.write("Warning: MGR webkey not defined\n")
#        s1 = "curl http://api.metagenomics.anl.gov/sequenceset/mgm%s-050-1/ -D /tmp/fetch-dump > %s.gz"  % ( a, a ) 
        s1 = "curl http://api.metagenomics.anl.gov/sequenceset/mgm%s-050-1/                     > %s.gz"  % ( a, a ) 
    else: 
        sys.stderr.write("Using MGR webkey %s\n" % key)
#        s1 = "curl 'http://api.metagenomics.anl.gov/sequenceset/mgm%s-050-1/?&auth=%s' -D /tmp/fetch-dump > %s.gz" % ( a, key, a )
        s1 = "curl 'http://api.metagenomics.anl.gov/sequenceset/mgm%s-050-1/?&auth=%s'  > %s.gz" % ( a, key, a )
#    jsonobject = os.popen(s).read()
#    print jsonobject
    sys.stderr.write("Executing %s\n" % s1) 
    os.popen(s1)
  
def retrieveSRRbyaccession(accession, rformat="fastq"):
    try:
        asperapath=os.environ["ASPERAPATH"]
    except KeyError:
        sys.exit("Error: environment variable ASPERAPATH must be defined and point to aspera to use SRR download")
    print  "Fetching SRR accession %s" % accession
    tla      = accession[0:3]
    stem     = accession[0:6]
    filename = accession+".sra"
    s = "%s/bin/ascp -l 300m -QT -i %s/etc/asperaweb_id_dsa.putty anonftp\@ftp-trace.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s ./%s" % ( asperapath, asperapath, tla, stem, accession, filename, filename)
    print "with %s" % (s)
    if int(os.popen("ascp 2>&1 |wc").read().split()[0]) > 10:
        os.system(s)
    else: 
        sys.exit("Sorry, can't find ascp.\t")
    print s 
    s = "fastq-dump --split-3 ./%s & " % (filename)
    print s
    os.system(s)
  
if __name__ == '__main__':
    usage  = "usage: fetch.py <accession number> [-f <format>]"
    parser = OptionParser(usage)
    parser.add_option("-f", "--format", dest="rformat", default=None, help="Data format (fasta, fastq, gbk) ")
    (opts, args) = parser.parse_args()
    try:
        key = os.environ["MGRKEY"]
    except KeyError:
        key = ""
  
    try :
        accession = args[0]
    except IndexError:
        parser.error("accession is a required parameter\n%s" % usage)
  
    if re.search("^(4......\..)$", accession):  # MGR accession
        retrieveMGRbyaccession(accession)
    elif re.search("^(MGR\d\d\d\d\d\d.\d)$", accession):  #  MGR accession 
        retrieveMGRbyaccession(accession, rformat=opts.rformat)
    elif re.search(".RR(......)", accession):   # SRR accession
        retrieveSRRbyaccession(accession)
    elif re.search("^(\w\w\w\w)$", accession):  # WGS accession
        retrieveWGSbyaccession(accession)
    elif re.search("^(N._\d\d\d\d\d\d)$", accession):  # REFSEQ accession
        retrieveGENBANKbyaccession(accession, rformat=opts.rformat)
    elif re.search("^(N._\d\d\d\d\d\d.\d)$", accession):  # REFSEQ accession
        retrieveGENBANKbyaccession(accession, rformat=opts.rformat)
    elif re.search("^(C.\d\d\d\d\d\d)$", accession):  #  Provisional REFSEQ
        retrieveGENBANKbyaccession(accession, rformat=opts.rformat)
    elif re.search("^(.U\d\d\d\d\d\d)$", accession):  #  Provisional REFSEQ
        retrieveGENBANKbyaccession(accession, rformat=opts.rformat)
    else: 
        print "Don't recognize acecssion %s" % accession
        sys.exit()

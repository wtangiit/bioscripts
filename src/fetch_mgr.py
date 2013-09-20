#!/usr/bin/env python
'''download MG-RAST input by accession number'''
import datetime
import json
import requests
import sys
import pprint
import time
from optparse import OptionParser

valid_file_id = ["050.1", "050.2"]  #fastq or fasta
               
def fetchByAccession(id, key=""):
    url = "http://api.metagenomics.anl.gov/1/download/mgm"+id
    if key != "":
        url += "?auth="+key
    
    r = requests.get(url, headers={}, allow_redirects=True)
    data = r.json()
    for item in data:
        if item.get("file_id") in valid_file_id:
            downloadurl =  item.get("url")
            if key != "":
                downloadurl += "&auth="+key
            
            filename = id+"."+item.get("file_name")
            print "downloading file from url %s to file %s..." % (downloadurl, filename)
            try:
                rget = requests.get(downloadurl, headers={}, allow_redirects=True, stream = True)
            except Exception as e:
                raise Exception(u'Unable to connect to MG-RAST API server %s: %s' %(downloadurl, e))
            if not (rget.ok):
                raise Exception(u'Unable to connect to MG-RAST API server %s: %s' %(downloadurl, rget.raise_for_status()))
            with open(filename, 'wb') as f:
                for chunk in rget.iter_content(chunk_size=8192): 
                    if chunk:
                        f.write(chunk)
                        f.flush()
            print "downloading done!"

if __name__ == "__main__":
    p = OptionParser()
    p.add_option("-k", dest = "key", type = "string", 
                    help = "MG-RAST webkey")
    
    p.add_option("-a", dest = "accession", type = "string", 
                    help = "MG-RAST accession id")
    
    (opts, args) = p.parse_args()
    
    if not opts.accession:
        print "please specify accession id, e.g. 4472152.3"
        p.print_help()
        exit()
        
    if opts.key:
        fetchByAccession(opts.accession, opts.key)
    else:
        fetchByAccession(opts.accession)

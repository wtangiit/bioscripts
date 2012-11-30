#!/usr/bin/env python

import sys, os
from optparse import OptionParser
import numpy as np

# import matplotlib.pyplot as plt

''' Distance matrix prototype -- loads abundance table (in COUNTS!) and returns all-agianst-all 
distance matrix.  Datasets in columns, taxa in rows; first column taxa headers, first row dataset headers.
Outputs an nxn matrix of floats to standard out'''

if __name__ == '__main__':
  usage  = "usage: %prog <input abundance table> [-d jsdbayes|euclidean ] \npurpose: computes all-against-all distance matrix using Jensen-Shannon divergence on Bayesian posteriors"
  parser = OptionParser(usage)
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose")
  parser.add_option("-d", "--distance",  dest="dist", default="euclidean", help="Distance metric: jsdbayes or euclidean")
  
  # parse options
  (opts, args) = parser.parse_args()
  if len(args) == 0 :
    parser.error("Input file is mandatory argument" )
   
  filename = args[0]
  if not (filename and os.path.isfile(filename) ):
    parser.error("Missing input file" )

  # Load datafile
  if opts.verbose: sys.stdout.write("Reading %s ... \n"%filename)
  datatable = np.loadtxt(filename, skiprows=1, dtype = "str", delimiter="\t")
  if opts.verbose:  print datatable
  firstcol = datatable[:,0]
  if opts.verbose: print firstcol
  datatable[np.nonzero(datatable == "")] = "0"       # normalize empty fields
  datatable = np.array(datatable[:,1:], dtype="float")  # chop off text column
  if opts.verbose:  print datatable
  if opts.verbose:  np.savetxt(sys.stdout, datatable, fmt="%d", delimiter="\t") 
  
  # evaluate distance matrix
  if opts.dist == "euclidean":
    A = datatable
    R = (A / np.sum(A.T,axis=1))       # divide by sum of each col
    R[np.nonzero(np.isnan(R)) ] = 0    # set nans to zero  (shouldn't be necessary)
#    plt.imshow(R, interpolation="nearest", aspect="auto")
#    plt.show()
    if opts.verbose:    np.savetxt(sys.stdout, R, fmt="%.4f", delimiter="\t")
    nsamples = len(A[0,:])
    distancemat= np.zeros((nsamples, nsamples)) 
    for i in range(0, nsamples  ) :
      for j in range(0,  i ) :
         d = R[:,i] - R[:,j]  
         D = np.sum(d*d)
         distancemat[i][j] = D
         distancemat[j][i] = D
    np.savetxt(sys.stdout, distancemat, fmt="%.4f", delimiter="\t") # display distancemat
  elif opts.dist == "jsdbayes":
    print "Not implemented"
  elif opts.dist == "bogot":
    print "Not implemented"
  else :
    print "Distance metric %s not recognized"%opts.dist


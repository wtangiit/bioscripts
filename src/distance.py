#!/usr/bin/env python

import sys, os
from optparse import OptionParser
import numpy as np
import scipy 
from scipy import stats

import matplotlib.pyplot as plt

''' Distance matrix prototype -- loads abundance table (in COUNTS!) and returns all-agianst-all 
distance matrix.  Datasets in columns, taxa in rows; first column taxa headers, first row dataset headers.
Outputs an nxn matrix of floats to standard out

euclidean -- euclidean distance, for testing
bogot     -- bogus t-distance.  A horrible approximation (gaussian approx to beta posteriors)
             WARNING: this is positive and symmetric but not metric

others are experimental and not warrantied.
'''

def normalize(y):
   k = y
   k[np.isinf( y ) ] =  0    # clean up 
   k[np.isnan( y ) ] =  0
   k = k / k.sum()
   return(k)

def kl2(a1, b1, a2, b2, x):    # KL divergence given x, log(posterior1) and log(posterior2)
   d = x[1] - x[0]
   p_lin = np.exp(p)
   p_lin[np.isinf( p_lin ) ] =  0    # clean up before sum
   p_lin[np.isnan( p_lin ) ] =  0
   p_lin = p_lin / p_lin.sum() / d   # Normalize
   p_lin[np.isinf( p_lin ) ] =  0    # clean up after divide
   p_lin[np.isnan( p_lin ) ] =  0
   print "pq, ", p
   print q
   ratio_pq = p - q
   ratio_pq[np.isinf( ratio_pq ) ] =  0
   a = (ratio_pq ) * p_lin * d
   a[np.isnan( a) ] =  0
   a[np.isinf( a) ] =  0
   return a.sum()

def kl(p, q, x):    # KL divergence given x, log(posterior1) and log(posterior2)
   d = x[1] - x[0]
   p_lin = np.exp(p)
   p_lin = p_lin / p_lin.sum() / d   # Normalize
   print "pq, ", p
   print q
   ratio_pq = p - q
   ratio_pq[np.isinf( ratio_pq ) ] =  0
   a = (ratio_pq ) * p_lin * d
   a[np.isinf( p_lin ) ] =  0  
   a[np.isnan( p_lin ) ] =  0
   a[np.isinf( p_lin ) ] =  0  
   a[np.isnan( p_lin ) ] =  0
   a[np.isnan( a) ] =  0
   a[np.isinf( a) ] =  0
   return a.sum()

def multid_kl(p, q, x):
   q[np.nonzero( q < 10**(-50)) ] =  10**(-50)
   assert q.all() > 0 
   d = x[1] - x[0]
   a = np.log(p / q ) * p * d
   a[np.isnan( a) ] =  0
#   print "a  ", a
#   print "kl  ", a.sum()
   return a.sum()

def multid_jsd(p, q):
   m = (p + q ) * 0.5
   b = multid_kl(p,m,x) * 0.5 + multid_kl(q, m, x) * 0.5
#   print "jsd ", b
   if b < 0 : b=0
   return b

def jsd(p, q,x):  
   d= x[1]- x[0]
   p_lin = np.exp(p)
   q_lin = np.exp(q)
   p_lin[np.isinf( p_lin ) ] =  0
   q_lin[np.isinf( q_lin ) ] =  0
   p_lin = p_lin / p_lin.sum() 
   q_lin = q_lin / q_lin.sum() 
   m = np.log((p_lin + q_lin)  * 0.5)
   if any(np.isnan(m)) : print "m nan"
   if any(np.isinf(m)) : print "m inf"
#   m[np.isnan( m ) ] =  -500
#   m[np.isinf( m ) ] =  -500
   b = kl(p,m,x) * 0.5 + kl(q, m,x) * 0.5
#   print "one %f, two %f  b %f"% (kl(p,m,x) * 0.5 , kl(q, m,x) * 0.5 , b) 
   if b < 0 : b=0
   return b

def showdistancemat(dist, headers):
   print "Dataset\t"+"\t".join(headers)
#   sys.stderr.write( "size %d, %d\n"%(len(dist[:,0]), len(dist[0,:] ) ) )
   for i in range(0,len(dist[:,0])) :
     print headers[i],
     for r in range(0, len(dist[i,:])) : 
        print "\t%.4f"%(dist[i,r]),

     print
  
def calculateposteriors(A, x, returnformat="linear"):
    posterior = np.zeros(( len(A[:,0]), nsamples, len(x) ))   # initialize
    T = np.sum(A.T,axis=1)                                    # sum each col
    for k in range(0,len(A[:,0])):   # iterate over features
      sys.stderr.write("k = "+str(k)+"\n" )
      for i in range(0, nsamples):   # construct posteriors, iterating over feature (k) and sample (i)  
#         print "beta of ",  datatable[k][i] + pseudo , T[i] - datatable[k][i] + pseudo , "T(i)  ", T[i]
         if returnformat=="linear":
           posterior[k][i] = scipy.stats.beta.pdf(x, datatable[k][i] + pseudo, T[i] - datatable[k][i] + pseudo  )
         elif returnformat=="logold":
           posterior[k][i] = scipy.stats.beta.logpdf(x, datatable[k][i] + pseudo, T[i] - datatable[k][i] + pseudo  )
         elif returnformat=="log":
           temp = (datatable[k][i] + pseudo -1 ) * np.log( x )  + ( T[i] - datatable[k][i] + pseudo - 1 ) *  np.log(1-x)
           temp2 = np.exp(temp)
           temp2 = temp2 / temp2.sum() 
            
           b  = np.log( temp2 )
#           print "sum: ", np.exp(b).sum()
#           print b
           posterior[k][i] = b 
    return posterior

 
# construct posteriors, iterating over feature (k) and sample (i) 
if __name__ == '__main__':
  usage  = "usage: %prog <input abundance table> [-d jsdbayes|euclidean|kl ] \npurpose: computes all-against-all distance matrix using Jensen-Shannon divergence on Bayesian posteriors"
  parser = OptionParser(usage)
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose")
  parser.add_option("-d", "--distance",  dest="dist", default="euclidean", help="Distance metric: jsdbayes or euclidean")
  parser.add_option("-p", "--pseudo",  dest="pseudo", default=0.1 , help="Regularization paramter ")

  # parse options
  (opts, args) = parser.parse_args()
  pseudo = float(opts.pseudo)
  if len(args) == 0 :
    parser.error("Input file is mandatory argument" )
   
  filename = args[0]
  if not (filename and os.path.isfile(filename) ):
    parser.error("Missing input file" )

  # Load datafile
  if opts.verbose: sys.stdout.write("Reading %s ... \n"%filename)
  fp = open(filename, "r")
  firstrow=fp.readline().rstrip()
  headers = firstrow.split("\t")[1:]
  fp.close()
  datatable = np.loadtxt(filename, skiprows=1, dtype = "str", delimiter="\t")
  if opts.verbose:  print datatable
  firstcol = datatable[:,0]
  if opts.verbose: print firstcol
  datatable[np.nonzero(datatable == "")] = "0"       # normalize empty fields
  datatable = np.array(datatable[:,1:], dtype="float")  # chop off text column
  if opts.verbose:  print datatable
  if opts.verbose:  np.savetxt(sys.stdout, datatable, fmt="%d", delimiter="\t") 
  
  nsamples = len(datatable[0,:])
  # evaluate distance matrix
  if opts.dist == "euclidean":
    A = datatable
    R = (A / np.sum(A.T,axis=1))       # divide by sum of each col
    R[np.nonzero(np.isnan(R)) ] = 0    # set nans to zero  (shouldn't be necessary)
    if opts.verbose:    np.savetxt(sys.stdout, R, fmt="%.4f", delimiter="\t")
    print "Input data matrix: (%d x %d) "%(len(A[:,0] ) , len(A[0,:]) )
    distancemat= np.zeros((nsamples, nsamples)) 
    for i in range(0, nsamples  ) :
      for j in range(0,  i ) :
         d = R[:,i] - R[:,j]  
         D = np.sum(d*d)
         distancemat[i][j] = D
         distancemat[j][i] = D
    showdistancemat(distancemat, headers)

  if opts.dist == "jsdnaive":
    A = datatable
    R = (A / np.sum(A.T,axis=1))       # divide by sum of each col
    R[np.nonzero(np.isnan(R)) ] = 0    # set nans to zero  (shouldn't be necessary)
    if opts.verbose:    np.savetxt(sys.stdout, R, fmt="%.4f", delimiter="\t")
    print "Input data matrix: (%d x %d) "%(len(A[:,0] ) , len(A[0,:]) )
    distancemat= np.zeros((nsamples, nsamples)) 
    
    for i in range(0, nsamples  ) : 
      for j in range(0,  i ) : 
         d=0
         m = 0.5 * ( R[:,i] + R[:,j] ) 
         for k in range(0, len( R[:,i] ) ) :
           if m[k] > 0 :
             print d, 0.5*  R[k,i] * (np.log(R[k,i]) - np.log(m[k])) + 0.5*  R[k,j] * (np.log(R[k,j]) - np.log(m[k]))
             d+= 0.5* R[k,i] * (np.log(R[k,i]) - np.log(m[k])) + 0.5*  R[k,j] * (np.log(R[k,j]) - np.log(m[k]))
         distancemat[j][i] = d/np.log(2)
    showdistancemat(distancemat, headers)

  elif opts.dist == "jsdbayes1":
    A = datatable
    T = np.sum(A.T,axis=1)       # sum each col
    x = np.arange(0,1 - .000005 ,.00001) + .000005
    logposterior    = calculateposteriors(A, x, returnformat="log")

    for k in [0,1]:                 # iterate over two features
      distancemat = np.zeros((nsamples, nsamples))
      for i in range(0, nsamples  ) :
         for j in range(0, nsamples ) :   
          distancemat[i][j] = np.sqrt(jsd(logposterior[k][i], logposterior[k][j], x ))
          distancemat[i][j] = np.sqrt(jsd(logposterior[k][i], logposterior[k][j], x ))
      showdistancemat(distancemat, headers)
  
  elif opts.dist == "kl1":
    A = datatable
    T = np.sum(A.T,axis=1)       # sum each col
    x = np.arange(0,1 - .000005 ,.00001) + .000005
    logposterior    = calculateposteriors(A, x, returnformat="log")
    for k in [0,1]:                 # iterate over two features
      distancemat = np.zeros((nsamples, nsamples))
      for i in range(0, nsamples  ) :
         for j in range(0, nsamples ) :  
          distancemat[i][j] = kl(logposterior[k][i], logposterior[k][j], x )
      showdistancemat(distancemat, headers)
  elif opts.dist == "kl2":
    A = datatable
    T = np.sum(A.T,axis=1)       # sum each col
    x = np.arange(0,1 - .000005 ,.00001) + .000005
    logposterior    = calculateposteriors(A, x, returnformat="log")
    for k in [0,1]:                 # iterate over two features
      distancemat = np.zeros((nsamples, nsamples))
      for i in range(0, nsamples  ) :
         for j in range(0, nsamples ) :  
          distancemat[i][j] = kl2(logposterior[k][i], logposterior[k][j], x )
      showdistancemat(distancemat, headers)


  elif opts.dist == "jsdbayes":
    A = datatable
    x = np.arange(0,1 - .005 ,.01) + .005
    posterior    = calculateposteriors(A, x, returnformat="linear")
    logposterior = calculateposteriors(A, x, returnformat="log")

    distancemat = np.zeros((nsamples, nsamples))
    for i in range(0, nsamples  ) :
         for j in range(0, i +1 ) :   
          distancemat[i][j] = np.sqrt(multid_jsd(posterior[:,i], posterior[:,j]) )
          distancemat[j][i] = np.sqrt(multid_jsd(posterior[:,i], posterior[:,j]) )
    showdistancemat(distancemat, headers)

  elif opts.dist == "bogot":
    A = datatable + pseudo
    T = np.sum(A.T,axis=1)             # sum each col after pseudocount regularization
    R = (A / np.sum(A.T,axis=1))       # divide by sum of each col
    R[np.nonzero(np.isnan(R)) ] = 0    # set nans to zero  (shouldn't be necessary)
    s2 = R * (1-R) / T

    distancemat = np.zeros((nsamples, nsamples))
    for i in range(0, nsamples  ) :
       for j in range(0, nsamples ) :  
        bogos2 = s2[:,i] + s2[:,j] 
        distancemat[i][j] = np.sqrt(( ( R[:,i] - R[:,j] ) ** 2  / bogos2 ).sum() )
    showdistancemat(distancemat, headers)
  else :
    print "Distance metric %s not recognized"%opts.dist


#!/usr/bin/env python

'''do linear fit for generated gene/rgene file

take trained gene or rgene as input, generate a new trained file gene.fit or rgene.fit where
transition probability is the linear function of gc content. In the mean time, generate the 
plots showing linear regression.

./linear_fit -i gene

'''
import os
from optparse import OptionParser
from numpy import arange,array,ones,linalg
from pylab import plot,show
import matplotlib.pyplot as plt

digit2nt = {0: 'A', 1:'C', 2:'G', 3:'T'}

linestyle1 = ['ob', 'or', 'ok', 'og']
linestyle2 = ['-b', '-r', '-k', '-g',]

def number2dimer(number):
    digit1 = number / 4
    digit2 = number % 4
    nt1 = digit2nt[digit1]
    nt2 = digit2nt[digit2]
    return nt1+nt2    

def list_fit(ylist, m, ij, weight_line_values=[]):
    xi = arange(0,45)
    A = array([ xi, ones(45)])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    from_dimer = number2dimer(ij)
    plt.title("M%d P(X|%s)" % (m, from_dimer))
    
    plots = []
    legend_labels = []
    
    parameter_list = []
    
    for k in range(4):
        y = ylist[k]
        to_nt = digit2nt[k]
        
        w = linalg.lstsq(A.T,y)[0] # obtaining the parameters
        parameter_list.append(w)
        
        line = w[0]*xi+w[1] # regression line
        plots.append(ax.plot(range(26, 71),line,linestyle2[k],range(26, 71),y, linestyle1[k]))
        legend_labels.append(to_nt)
        
        print "plotting line M%s:P(%s|%s), slope=%.4f, inception=%.4f" % (m, to_nt, from_dimer, w[0], w[1])
        
    if len(weight_line_values) > 0: #plot weight line
        ax2 = ax.twinx()
        plots.append(ax2.plot(range(26, 71), weight_line_values, '-xm', markeredgewidth=1))
        ax2.set_ylim(0, 400000)
        ax2.set_ylabel('triplet count')
        
        
    ax.legend((plots[0][1], plots[1][1], plots[2][1], plots[3][1], ), ('A', 'C', 'G', 'T'), loc=0)
            
    ax.set_xlabel('GC content')
    ax.set_ylabel('probability')
    ax.grid(True)
    
    savefile = "M%d_%s.png" % (m, from_dimer)
    plt.savefig(savefile)
    
    return parameter_list
            
def parse_file(filename):
    infile = open(filename, "r")

    data_lists = [[[[] for k in range(4)] for ij in range(16)] for m in range(6)]# 6x16x4x[45]
        
    for g in range(45):
        line = infile.readline().strip('\n')
        for m in range(6):
            for ij in range(16):
                line = infile.readline().strip('\n')
                numbers = line.split()
                for k in range(4):
                    prob = numbers[k]
                    data_lists[m][ij][k].append(prob)
    return data_lists

def gen_fitted_file(linear_parameters):
    '''generated a new trained gene transition model with fitted linear functions'''
    
    fit_data_lists = [[[[] for ij in range(16)] for m in range(6)] for g in range(45)] #45x6x16x4
    
    for m in range(6):
        for ij in range(16):
            for k in range(4):
                w = linear_parameters[m][ij][k]
                slope = w[0]
                intercept = w[1]
                for g in range(45):
                    if g < 5:
                        xi = 5
                    else:
                        xi = g  
                    prob = slope * xi + intercept
                    if prob < 0.0001:
                        prob = 0.0001
                        
                    fit_data_lists[g][m][ij].append(prob)
    
    outfile = open("gene.fit", "w")
    
    for g in range(45):
        gc = g+26
        line = "%s\n" % gc
        outfile.write(line)
        for m in range(6):
            for ij in range(16):
                line = ""
                for k in range(4):
                    prob = fit_data_lists[g][m][ij][k]
                    line += "%.4f\t" % prob
                line.strip('\t')
                line += '\n'
                outfile.write(line)
    outfile.close()    
    
if __name__ == '__main__':
    usage  = "usage: %prog -i <input training data>"
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",  dest="input", type = "string", default=None, help="<input training data>")
    parser.add_option("-w", "--weight",  dest="weight", type = "string", default=None, help="<input weighting data>")
    
    (opts, args) = parser.parse_args()
    
    data_lists = parse_file(opts.input)  #datalists array of 6x16x4x45
    
    if opts.weight:
        weight_lists = parse_file(opts.weight)
            
    linear_parameter = [[[]for ij in range(16)] for m in range(6)]
    
    for m in range(6):
        for ij in range(16):
            ylist = data_lists[m][ij]    #ylist  4x45
            
            weight_line = []
            
            if opts.weight:
                for g in range(45):
                    weight_ct = weight_lists[m][ij][0][g]
                    weight_line.append(weight_ct)
                        
            parameters = list_fit(ylist, m, ij, weight_line)
            
            
            linear_parameter[m][ij] = parameters
                
    gen_fitted_file(linear_parameter)
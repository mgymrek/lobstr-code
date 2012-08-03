#!/usr/bin/python

def usage():
    print """
./plotAllelotypes.py <genotypes.tab plot debug file> <outputprefix>
"""
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

try:
    f = open(sys.argv[1],"r")
    prefix = sys.argv[2]
except:
    usage()
    sys.exit(1)

# get each locus
line = f.readline()
locusstring = ""
while line != "":
    while "Posterior" not in line and line != "":
        if "ProcessLocus" in line:
            locusstring = ".".join(line.strip().split()[1:])
        line = f.readline()
    if line == "": break
    # fill in numpy grid of posterior
    minal = int(line.strip().split()[6])
    maxal = int(line.strip().split()[7])
    datarange = maxal-minal+1
    data = [0.0 for i in xrange(int(datarange*datarange))]
    data = np.array(data)
    data.shape = (datarange,datarange)
    while "Posterior" in line:
        items = line.split()
        a1 = int(items[3])
        a2 = int(items[4])
        score = float(items[5])
        data[a1-minal, a2-minal] = score
        data[a2-minal, a1-minal] = score
        line = f.readline()
    print locusstring
    plt.clf()
    fig = plt.figure()
    plt.pcolormesh(data,cmap="Blues")
    plt.colorbar()
    plt.savefig("%s.%s.png"%(prefix,locusstring))
    if "ProcessLocus" in line:
        locusstring = ".".join(line.strip().split()[1:])
    line = f.readline()
    if line.strip() == "": break

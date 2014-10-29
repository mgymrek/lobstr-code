#!/usr/bin/env python
"""
This script is run on output of the
allelotype tool run with --command train. It generates an HTML file providing
visualization of PCR stutter noise at STRs.
"""
import matplotlib
matplotlib.use('Agg')
import argparse
import StringIO
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import urllib, base64
from scipy.stats.mstats import mquantiles
pd.options.mode.chained_assignment = None

NOISE_MODEL_PREFIX = None
OUT_FILE = None
MAX_PERIOD = 6

class StepModel(object):
    """
    Contains information about step size distributions
    """
    def __init__(self, filename):
        self.LoadStepModel(filename)

    def LoadStepModel(self, filename):
        self.NonUnitStepByPeriod = {}
        self.StepSizeByPeriod = {}
        f = open(filename, "r")
        for i in range(MAX_PERIOD):
            line = f.readline()
            self.NonUnitStepByPeriod[i+1] = float(line.strip())
        self.ProbIncrease = float(f.readline().split("=")[1].strip())
        for i in range(MAX_PERIOD):
            line = f.readline()
            items = map(float, line.strip().split()[1:])
            self.StepSizeByPeriod[i+1] = items
        f.close()

    def PlotHistogram(self, period):
        """ Plot histogram of error sizes"""
        colors = ["gray","red","gold","blue","green","purple"]
        fig = plt.figure()
        fig.set_size_inches((6,3))
        ax = fig.add_subplot(111)
        xran = range(-18,19)
        ax.bar(xran, self.StepSizeByPeriod[period], align="center", color=colors[period-1], alpha=0.3)
        ax.bar([xran[i] for i in xran if xran[i]%period==0], [self.StepSizeByPeriod[period][i] for i in xran if xran[i]%period==0], align="center", color=colors[period-1])
        ax.set_xlabel("Step size", size=15)
        ax.set_ylabel("Frequency", size=15)
        ax.set_xlim(left=-3*MAX_PERIOD, right=3*MAX_PERIOD)
        ax.set_yticklabels(ax.get_yticks(), size=15)
        ax.set_xticks(xran)
        ax.set_xticklabels(xran, size=12, rotation=90)
        imgdata = StringIO.StringIO()
        fig.savefig(imgdata, format="png")
        imgdata.seek(0)
        return urllib.quote(base64.b64encode(imgdata.buf))

class StutterProblem(object):
    """ Contains info on reads used to build stutter model """
    def __init__(self, filename):
        self.LoadProblem(filename)

    def LoadProblem(self, filename):
        sp = pd.read_csv(filename, sep=" ", names=["stutter","period","length","gc","score"])
        sp = sp[sp["period"].apply(lambda x: type(x)==str)] # get rid of empty last line
        sp.loc[:,"period"] = sp["period"].apply(lambda x: int(x.split(":")[1]))
        sp.loc[:,"length"] = sp["length"].apply(lambda x: int(x.split(":")[1]))
        sp.loc[:,"gc"] = sp["gc"].apply(lambda x: float(x.split(":")[1]))
        sp.loc[:,"score"] = sp["score"].apply(lambda x: float(x.split(":")[1]))
        self.sp = sp

    def GetAverageStutterProb(self):
        return self.sp[self.sp["stutter"]==1].shape[0]*1.0/self.sp.shape[0]

    def StutterByPeriod(self, period):
        tmp = self.sp[self.sp.period==period]
        return tmp[tmp.stutter==1].shape[0]*1.0/tmp.shape[0]

    def PlotVariable(self, variable, bins):
        colors = ["gray","red","gold","blue","green","purple"]
        fig = plt.figure()
        fig.set_size_inches((6,3))
        ax = fig.add_subplot(111)
        for period in range(1, MAX_PERIOD+1):
            sub = self.sp[self.sp.period==period]
            # get deciles of the variable of interest
            probs = []
            for i in range(len(bins)-1):
                l = bins[i]
                u = bins[i+1]
                tmp = sub[(sub[variable]>=l) & (sub[variable]<u)]
                probs.append(np.mean(tmp.stutter+1)*0.5)
            ax.plot(bins[:-1], probs, color=colors[period-1], label="Period %s"%period)
        ax.set_xlabel(variable, size=15)
        ax.set_ylabel("Stutter probability", size=15)
        ax.set_xticklabels(ax.get_xticks())
        ax.set_yticklabels(ax.get_yticks())
        ax.legend(loc="upper left")
        imgdata = StringIO.StringIO()
        fig.savefig(imgdata, format="png")
        imgdata.seek(0)
        return urllib.quote(base64.b64encode(imgdata.buf))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--noise_model", help="Prefix to noise model output from allelotype.", type=str, required=True)
    parser.add_argument("--out", help="Name of output html file", type=str, required=True)
    args = parser.parse_args()
    NOISE_MODEL_PREFIX = args.noise_model
    OUT_FILE = args.out

    # Read in noise model files
    stepmodel = StepModel(NOISE_MODEL_PREFIX + ".stepmodel")
    stutterproblem = StutterProblem(NOISE_MODEL_PREFIX + ".stutterproblem")

    # Set up HTML file
    f = open(OUT_FILE, "w")
    f.write("<head><title>lobSTR PCR stutter analysis - %s</title></head>\n"%NOISE_MODEL_PREFIX)
    f.write("<body>\n")

    # Write stats from noise model
    f.write("<h1>Stutter params</h1>\n")
    f.write("<h2>Stutter probability by period</h2>\n")
    for i in range(MAX_PERIOD): f.write("Period %s: %s<br>\n"%(i+1, stutterproblem.StutterByPeriod(i+1)))
    f.write("<br>")
    f.write("Average stutter rate: %s<br>\n"%(stutterproblem.GetAverageStutterProb()))
    f.write("Percent of stutter errors that increase the repeat number: %s<br>\n"%(stepmodel.ProbIncrease))

    # Stutter err distributions
    f.write("<h1>Error size distributions</h1>\n")
    for i in range(MAX_PERIOD):
        f.write("<h2>Period %s</h2>\n"%(i+1))
        pltdata = stepmodel.PlotHistogram(i+1)
        uri = "data:image/png;base64," + pltdata
        f.write("<img src=\"%s\"/>\n"%uri)

    # Stutter by feature (track length, GC, purity, by period)
    f.write("<h1>Sequence features</h1>\n")
    featureBins = {"length": np.arange(10, 60, 10),
                   "score": [0.5, 0.8, 1],
                   "gc": np.arange(0.3,0.7,0.1)}
    for feature in ["length","gc","score"]:
        f.write("<h2>%s</h2>\n"%feature)
        pltdata = stutterproblem.PlotVariable(feature, featureBins[feature])
        uri = "data:image/png;base64," + pltdata
        f.write("<img src=\"%s\"/>\n"%uri)


    f.write("</body>\n")
    f.close()

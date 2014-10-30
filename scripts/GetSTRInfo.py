#!/usr/bin/env python
"""python GetSTRInfo.py <STR table> <ref.fa>

Prints output to stdout

For each locus, get:
TRF score
GC content
Entropy of 50bp up/downstream
"""

import sys
import math
from Bio.Seq import Seq
from Bio import SeqIO

FLANKLEN = 50

def loadGenome(genomeFile):
    """ load human genome"""
    handle = open(genomeFile,"rU")
    genome = SeqIO.to_dict(SeqIO.parse(handle,"fasta"))
    return genome

def ExtractFlanks(genome, chrom, start, end):
    """ extract +/- FLANKLEN from STR locus """
    leftFlank = str(genome[chrom][max([start-FLANKLEN,0]):start].seq).upper()
    rightFlank = str(genome[chrom][end:min([len(genome[chrom]),end+FLANKLEN])].seq).upper()
    return leftFlank+rightFlank

def GetGC(flanks):
    """ get gc content of nucleotide string """
    gc = 0
    total = 0
    for i in flanks:
        if i != "N":
            total += 1
            if i == "C" or i == "G": gc += 1
    return gc*1.0/total

def GetEntropy(flanks):
    """ get sequence entropy of nucleotide string """
    countA = 0
    countT = 0
    countG = 0
    countC = 0
    for i in flanks:
        if i == "A":
            countA += 1
        elif i == "T":
            countT += 1
        elif i == "C":
            countC += 1
        elif i == "G":
            countG += 1
        else: pass
    total = countA+countT+countG+countC
    fractions = [item*1.0/total for item in [countA,countT,countG,countC]]
    entropy = sum([-1.0*item*math.log(item,2) for item in fractions if item != 0])
    return entropy

def ProcessEachLocus(strtablefile,genome):
    """ get entropy, gc, and score for each locus """
    f = open(strtablefile,"r")
    line = f.readline()
    while line != "":
        items = line.strip().split("\t")
        chrom = items[0]
        start = int(items[1])
        end = int(items[2])
        period = int(items[3])
        length = float(items[4])
        maxscore = int(period*length)*2
        motif = items[-1].strip()
        if maxscore > 0:
            score = float(items[8])/maxscore
            try:
                flanks = ExtractFlanks(genome,chrom,start,end)
                gc = GetGC(flanks)
                entropy = GetEntropy(flanks)
                print "\t".join(map(str,[chrom,start,end,score,gc,entropy]))
            except: pass
        line = f.readline()
    f.close()

def main():
    try:
        strtablefile = sys.argv[1]
        if strtablefile == "-": strtablefile = "/dev/stdin"
        genomeFile = sys.argv[2]
        if genomeFile == "-": genomeFile = "/dev/stdin"
    except:
        print __doc__
        sys.exit(1)
    sys.stderr.write("loading genome...\n")
    genome = loadGenome(genomeFile)

    sys.stderr.write("processing each locus...\n")
    print "\t".join(["chrom","start","end","score","GC","entropy"])
    ProcessEachLocus(strtablefile, genome)

main()

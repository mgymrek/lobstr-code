def usage():
    print """
python lobstr_index.py [OPTIONS] --str <str_file> --ref <ref_file> --out_dir <output directory>

OPTIONS:
--str: bed file containing str information
--ref: reference genome in fasta format
--out_dir: path to write results to
-h, --help: print this usage screen
-v: verbose
--extend: the length of flanking region to include on either side of the STR. (default 1000). Note, if you change this paramter, YOU MUST USE THE SAME --extend OPTION WHEN YOU CALL lobSTR

In $out_dir, creates:
Creates BWT references for each repeat unit
Also outputs a file with the sizes of chromosomes and a
table with all STR repeat units

Requires lobSTRIndex to be installed in the $PATH.
"""

import os
import sys
import getopt
from Bio.Seq import Seq
from Bio import SeqIO

try:
    opts, args = getopt.getopt(sys.argv[1:], "hv", ["help","str=", "ref=", "out_dir=", "extend="])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)
args = [item[0] for item in opts]

if (not("--str" in args and "--ref" in args and "--out_dir" in args)):
    usage()
    sys.exit(2)
# initialize variables
strfile = ""
reffile = ""
outdir = ""
extend = 1000
pad = 50
verbose = False

for o,a in opts:
    if o == "-v": verbose = True
    if o == "--str": strfile = a
    if o == "--ref": reffile = a
    if o == "--out_dir": outdir = a
    if o == "--extend": extend = int(a)
    if o == "-h" or o == "--help":
        usage()
        sys.exit(0)

###########################
# methods
nucToNumber={"A":0,"C":1,"G":2,"T":3}
def loadGenome(genomeFile):
    """ load human genome"""
    handle = open(genomeFile,"rU")
    genome = SeqIO.to_dict(SeqIO.parse(handle,"fasta"))
    return genome

def getCanonicalMS(repseq):   
    size = len(repseq)
    canonical = repseq
    for i in range(size):
        newseq = repseq[size-i:]+repseq[0:size-i]
        for j in range(size):
            if nucToNumber[newseq[j]] < nucToNumber[canonical[j]]:
                canonical = newseq
            elif nucToNumber[newseq[j]] > nucToNumber[canonical[j]]:
                break
    return canonical

def compareString(seq1,seq2):
    size = len(seq1)
    for i in range(size):
        if nucToNumber[seq1[i]] < nucToNumber[seq2[i]]: 
            return seq1
        if nucToNumber[seq1[i]] > nucToNumber[seq2[i]]:
            return seq2
    return seq1
        

def reverseComplement(seq):
    newseq = ""
    size = len(seq)
    for i in range(len(seq)):
        char = seq[len(seq)-i-1]
        if char == "A": newseq += "T"
        if char == "G": newseq +="C"
        if char == "C": newseq +="G"
        if char == "T": newseq += "A"
    return newseq

def PadFlank(seq, n):
    return "N"*n+seq+"N"*n

def processTRF(strfile, outdir, genome):
    repToFile = {} # map repseq -> file
    sevenmerdict = {}
    g = open(outdir+"/lobSTR_strdict.txt","w")
    # write chrom sizes to strdict
    for chrom in genome:
        print chrom, str(len(genome[chrom]))
        g.write("\t".join(["REF",chrom, str(len(genome[chrom]))])+"\n")

    f = open(strfile, "r")
    allfasta = open(outdir+"/lobSTR_ref.fa", "w")
    line = f.readline()
    ident = 0
    while line !="":
        items = line.strip().split("\t")
        chrom, start, end, copynum, repseq = items[0], int(items[1]), int(items[2]), items[4], items[14]
        try:
            name = items[15]
            if len(name) > 15: name = name[0:15] # truncate if too long
        except: name="."
        # extract flanking regions
        if "$" not in chrom:
            try:
                leftFlank = str(genome[chrom][max(start-extend,0):end].seq).upper()
                rightFlank = str(genome[chrom][start:min(end+extend,len(genome[chrom]))].seq).upper()
                strregion = str(genome[chrom][max(start-extend,0):min(end+extend,len(genome[chrom]))].seq).upper()
                repseq = getCanonicalMS(repseq)
                revrepseq = getCanonicalMS(reverseComplement(repseq))
                repseq = compareString(repseq, revrepseq)
            except:
                repseq = ""
        
        if len(repseq) <= 6 and len(repseq) >= 1 and (start-extend) > 0 and "$" not in chrom and len(strregion) > 0:
            # write fasta entries
            lident = ">"+"$".join(map(str,[ident,"L",chrom,start-extend,end,repseq, copynum, name]))
            rident = ">"+"$".join(map(str,[ident,"R",chrom,start,end+extend,repseq, copynum, name]))
            strident = ">"+"$".join(map(str,[ident,chrom,start-extend,end+extend,repseq,copynum,name]))
            try:
                repToFile[repseq].write(rident+"\n")
                repToFile[repseq].write(PadFlank(rightFlank,pad)+"\n")
                repToFile[repseq].write(lident+"\n")
                repToFile[repseq].write(PadFlank(leftFlank,pad)+"\n")
            except:
                repToFile[repseq] = open(outdir+"/lobSTR_%s.fa"%repseq,"w")
                repToFile[repseq].write(rident+"\n")
                repToFile[repseq].write(PadFlank(rightFlank,pad)+"\n")
                repToFile[repseq].write(lident+"\n")
                repToFile[repseq].write(PadFlank(leftFlank,pad)+"\n")
                # write ms dict entry
                g.write(repseq + "\n")
            # write whole region to master fasta reference
            allfasta.write(strident+"\n")
            allfasta.write(strregion+"\n")
            ident = ident + 1
        elif len(repseq) == 7 and "$" not in chrom:
             # write fasta entries
            lident = ">"+"$".join(map(str,[ident,"L",chrom,start-extend,end,repseq, copynum, name]))
            rident = ">"+"$".join(map(str,[ident,"R",chrom,start,end+extend,repseq, copynum, name]))
            if repseq not in sevenmerdict.keys():
                sevenmerdict[repseq] =(rident+"\n")
                sevenmerdict[repseq]+=(rightFlank+"\n")
                sevenmerdict[repseq]+=(lident+"\n")
                sevenmerdict[repseq]+=(leftFlank+"\n")
                # write ms dict entry
                g.write(repseq + "\n")
            else:
                sevenmerdict[repseq]+=(rident+"\n")
                sevenmerdict[repseq]+=(rightFlank+"\n")
                sevenmerdict[repseq]+=(lident+"\n")
                sevenmerdict[repseq]+=(leftFlank+"\n")
            ident = ident + 1
        line = f.readline()
    # close all the files
    f.close()
    g.close()
    for item in repToFile: repToFile[item].close()
    for item in sevenmerdict:
        f = open(outdir+"/lobSTR_%s.fa"%item,"w")
        f.write(sevenmerdict[item])
        f.close()
    allfasta.close()
    print "Processed %s records"%ident
    return repToFile, sevenmerdict

def bwaIndex(repToFile, sevenmerdict):
    for rep in repToFile:
        filename = repToFile[rep].name
        cmd = "lobSTRIndex index -a is %s"%filename
        os.system(cmd)
    for rep in sevenmerdict:
        filename = outdir+"/lobSTR_%s.fa"%rep
        cmd = "lobSTRIndex index -a is %s"%filename
        os.system(cmd)

###########################
def main():
    # load genome
    if (verbose): print "Loading genome..."
    genome = loadGenome(reffile)

    # make individual fasta files for each repeat unit
    # and write ms info to a file
    if (verbose): print "Processing STR table..."
    repToFile, sevenmerdict = processTRF(strfile, outdir, genome)
    
    # call bwa index on each of them
    if (verbose): print "Building bwa index..."
    bwaIndex(repToFile, sevenmerdict)

main()



DESCRIPTION = """
Create lobSTR index

In $out_dir, creates:
A bed file with merged reference regions
BWT index
A file with chromosome sizes
Table mapping reference regions to which STRs they contain

***Requires lobSTRIndex to be installed in the $PATH.***
"""

import argparse
import os
import pyfasta
import sys

STRREFFILE = None
REFFASTA = None
OUTDIR = None
EXTEND = 1000
VERBOSE = False
PAD = 50

###########################
# methods
nucToNumber={"A":0,"C":1,"G":2,"T":3}

def PROGRESS(msg):
    sys.stderr.write(msg.strip() + "\n")

def getCanonicalMS(repseq):
    """ Get canonical STR sequence """
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
    """ Compare two strings alphabetically """
    size = len(seq1)
    for i in range(size):
        if nucToNumber[seq1[i]] < nucToNumber[seq2[i]]: 
            return seq1
        if nucToNumber[seq1[i]] > nucToNumber[seq2[i]]:
            return seq2
    return seq1
        
def reverseComplement(seq):
    """ Get the reverse complement of a nucleotide string """
    newseq = ""
    size = len(seq)
    for i in range(len(seq)):
        char = seq[len(seq)-i-1]
        if char == "A": newseq += "T"
        if char == "G": newseq +="C"
        if char == "C": newseq +="G"
        if char == "T": newseq += "A"
    return newseq

def WriteChromSizes(genome, outfile):
    f = open(outfilel, "w")
    for chrom in genome.keys():
        f.write("\t".join([chrom, len(genome[chrom])])+"\n")
    f.close()

def PadFlank(seq, n):
    return "N"*n+seq+"N"*n

def GenerateMergedReferenceBed(str_ref_file, merged_str_file):
    """ Generated merged reference region """
    # TODO

def GetRefFasta(merged_str_file, str_ref_fasta, str_map_file):
    """ Get reference fasta and map file from the bed file """
    # TODO

def bwaIndex(str_ref_fasta):
    cmd = "lobSTRIndex index -a is %s"str_ref_fasta
    os.system(cmd)

##########################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("--str", help="Bed file containing STR information.", required=True, type=str)
    parser.add_argument("--ref", help="Reference genome in fasta format", required=True, type=str)
    parser.add_argument("--out_dir", help="Path to write results to", required=True, type=str)
    parser.add_argument("--extend", help="Length of flanking region to include on either side of the STR. (default 1000)", required=False, type=int)
    parser.add_argument("--verbose", help="Print out useful messages", required=False, action="store_true")

    args = parser.parse_args()
    STRREFFILE = args.str
    REFFILE = args.ref
    OUTDIR = args.out_dir
    if args.verbose: VERBOSE = True

    if VERBOSE: PROGRESS("Loading gnome...")
    genome = pyfasta.Fasta(REFFILE)
    WriteChromSizes(genome, os.path.join(OUTDIR, "lobSTR_chromsizes.tab"))

    if VERBOSE: PROGRESS("Processing STR table...")
    merged_str_file = os.path.join(OUTDIR, "lobSTR_mergedref.bed")
    str_ref_fasta = os.path.join(OUTDIR, "lobSTR_ref.fasta")
    str_map_file = os.path.join(OUTDIR, "lobSTR_ref_map.tab")
    GenerateMergedReferenceBed(STRREFFILE, merged_str_file)
    GetRefFasta(merged_str_file, str_ref_fasta, str_map_file)

    if VERBOSE: PROGRESS("Building BWA index...")
    bwaIndex(str_ref_fasta)

#!/usr/bin/env python
"""Create lobSTR index

In $out_dir, creates:
A bed file with merged reference regions
BWT index
A file with chromosome sizes
Table mapping reference regions to which STRs they contain

***Requires lobSTRIndex and bedtools to be installed in the $PATH.***

Note: this was tested with bedtools v2.22.1-17-gd6547b3. Older versions
may not be compatible with this script.
"""

import argparse
import os
import pyfasta
import shutil
from subprocess import Popen, PIPE, STDOUT
import sys
import tempfile

STRREFFILE = None
REFFASTA = None
OUTDIR = None
EXTEND = 1000
VERBOSE = False
PAD = 50
DEBUG = False

###########################
# methods
nucToNumber={"A":0,"C":1,"G":2,"T":3}

def PROGRESS(msg):
    sys.stderr.write(msg.strip() + "\n")

def ERROR(msg):
    sys.stderr.write(msg.strip() + "\n")
    sys.exit(1)

def RunCommand(cmd):
    if DEBUG:
        PROGRESS("CMD: %s"%cmd)
        return
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, \
                  stderr=STDOUT, close_fds=True)
    ex = p.wait()
    if ex != 0:
        stdout, stderr = "", ""
        if p.stdout is not None: stdout = p.stdout.read()
        if p.stderr is not None: stderr = p.stderr.read()
        PROGRESS("ERROR: command '%s' failed.\n\nSTDOUT:%s\nSTDERR:%s"%(cmd, stdout, stderr))
        sys.exit(1)

def IsExec(binname):
    """ Check if a file is executable """
    return os.path.isfile(binname) and os.access(binname, os.X_OK)

def CheckBin(binname):
    """ Check if a binary exists on the PATH """
    for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        exe_file = os.path.join(path, binname)
        if IsExec(exe_file):
            return True
    return False


def GetRepseq(repseq):
    """ Get canonical STR sequence, considering both strands """
    repseq_f = getCanonicalMS(repseq)
    repseq_r = getCanonicalMS(reverseComplement(repseq))
    repseq = compareString(repseq_f, repseq_r)
    return repseq

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
        if char == "G": newseq += "C"
        if char == "C": newseq += "G"
        if char == "T": newseq += "A"
    return newseq

def WriteChromSizes(genome, outfile):
    f = open(outfile, "w")
    for chrom in genome.keys():
        f.write("\t".join(map(str,[chrom.split()[0], len(genome[chrom])]))+"\n")
    f.close()

def PadFlank(seq, n):
    return "N"*n+seq+"N"*n

def GenerateMergedReferenceBed(str_ref_file, annotated_ref):
    """ Generated merged reference region """
    # Get temp file names
    tmpdir = tempfile.mkdtemp()
    sorted_ref_withflank = os.path.join(tmpdir, "lobSTR_ref_plusflank_sorted.bed")
    merged_ref = os.path.join(tmpdir, "lobSTR_ref_merged.bed")
    # Generate reference file
    cmd_sort = """cat %s | awk '{print $1 "\\t" $2-%s "\\t" $3+%s "\\t" $0}' | cut -f 4 --complement | \
awk '($2>0)' | cut -f 1-5,17 | sortBed -i stdin > %s """%(str_ref_file, EXTEND, EXTEND, sorted_ref_withflank)
    RunCommand(cmd_sort)
    cmd_merge = """mergeBed -i %s > %s"""%(sorted_ref_withflank, merged_ref)
    RunCommand(cmd_merge)
    cmd_annot = """intersectBed -a %s -b %s -wa -wb | awk '{print $1 "\\t" $2 "\\t" $3 "\\t" $7"_"$8"_"$9";"}' | \
bedtools groupby -g 1,2,3 -c 4 -o concat > %s"""%(merged_ref, sorted_ref_withflank, annotated_ref)
    RunCommand(cmd_annot)
    # Clean up
    shutil.rmtree(tmpdir)

def GetRefFasta(genome, refkeys, merged_str_file, str_ref_fasta, str_map_file):
    """ Get reference fasta and map file from the bed file """
    f_merged = open(merged_str_file, "r")
    f_fa = open(str_ref_fasta, "w")
    f_map = open(str_map_file, "w")
    # Write fasta and map entry for each reference chunk
    line = f_merged.readline()
    refnum = 0
    while line != "":
        chrom, start, end, annot = line.strip().split("\t")
        try:
            refchrom = refkeys[chrom]
        except:
            ERROR("Chromosome %s not in reference fasta\n"%chrom)
        start = int(start)
        end = int(end)
        if end >= len(genome[refchrom]): end = len(genome[refchrom])-1
        refseq = PadFlank(genome[refchrom][start:end], PAD).upper()
        annotations = annot.split(";")[:-1]
        motifs = [GetRepseq(item.split("_")[2]) for item in annotations]
        annotations = [annotations[i] + "_" + motifs[i] for i in range(len(annotations))]
        # Write to fasta file
        f_fa.write(">%s$%s$%s$%s\n"%(refnum, chrom, start, end))
        f_fa.write("%s\n"%refseq)
        # Write to map file
        f_map.write("\t".join(map(str,[refnum, ";".join(set(motifs)), ";".join(annotations)]))+"\n")
        # Get next record
        line = f_merged.readline()
        refnum += 1
    f_merged.close()
    f_fa.close()
    f_map.close()

def bwaIndex(str_ref_fasta):
    cmd = "lobSTRIndex index -a is %s"%str_ref_fasta
    RunCommand(cmd)

##########################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--str", help="Bed file containing STR information.", required=True, type=str)
    parser.add_argument("--ref", help="Reference genome in fasta format", required=True, type=str)
    parser.add_argument("--out_dir", help="Path to write results to", required=True, type=str)
    parser.add_argument("--extend", help="Length of flanking region to include on either side of the STR. (default 1000)", required=False, type=int, default=1000)
    parser.add_argument("--verbose", help="Print out useful messages", required=False, action="store_true")
    parser.add_argument("--debug", help="Don't run commands, just pring them", required=False, action="store_true")

    args = parser.parse_args()
    STRREFFILE = args.str
    REFFILE = args.ref
    OUTDIR = args.out_dir
    EXTEND = args.extend
    if not os.path.exists(OUTDIR):
        try:
            os.mkdir(OUTDIR)
        except OSError:
            ERROR("Could not create outdirectory %s"%OUTDIR)
    if args.verbose: VERBOSE = True
    if args.debug: DEBUG = True

    if not CheckBin("lobSTRIndex"): ERROR("Could not find lobSTRIndex on the $PATH. Please install lobSTR")
    if not CheckBin("sortBed"): ERROR("Could not find sortBed on the $PATH. Please install bedtools")
    if not CheckBin("mergeBed"): ERROR("Could not find mergeBed on the $PATH. Please install bedtools")
    if not CheckBin("intersectBed"): ERROR("Could not find intersectBed on the $PATH. Please install bedtools")
    if not CheckBin("bedtools"): ERROR("Could not find bedtools on the $PATH. Please install bedtools")

    if VERBOSE: PROGRESS("Loading genome...")
    genome = pyfasta.Fasta(REFFILE)
    refkeys = dict([(key.split()[0], key) for key in genome.keys()])
    if not DEBUG: WriteChromSizes(genome, os.path.join(OUTDIR, "lobSTR_chromsizes.tab"))

    if VERBOSE: PROGRESS("Processing STR table...")
    merged_str_file = os.path.join(OUTDIR, "lobSTR_mergedref.bed")
    str_ref_fasta = os.path.join(OUTDIR, "lobSTR_ref.fasta")
    str_map_file = os.path.join(OUTDIR, "lobSTR_ref_map.tab")
    GenerateMergedReferenceBed(STRREFFILE, merged_str_file)
    if not DEBUG: GetRefFasta(genome, refkeys, merged_str_file, str_ref_fasta, str_map_file)

    if VERBOSE: PROGRESS("Building BWA index...")
    bwaIndex(str_ref_fasta)

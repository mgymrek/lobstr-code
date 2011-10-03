def usage():
    print """
python strangers.py [OPTIONS] -f <file1[,file2,...]> -t <table filename> -g <genome.fa> -o <output prefix> [-d <database file> | -N]

-f,--files     file or comma-separated list of files containing reads in fasta or fastq format
-t,--table     file containing table of markers to test for. See README for table format
-g,--genome    fasta file containing the genome (e.g. hg18.fa)
-o,--out       prefix for out put files. will output:
                      <prefix>.aligned.tab
                      <prefix>.genotypes.tab
                      <prefix>.filename.ids
                      <prefix>.filename
                      <prefix>.haplotype
                      
-d, --database allele database with alleles for the anonymized haplotype. See readme. Either -d or -N must be specified.
-N             Mask revealing reads instead of anonymizing them.

Options:
-h,--help                  display this help screen
-v,--verbose               print out useful progress messages
-q,--fastq                 reads are in fastq format (default: fasta)
-p,--threads <int>         number of threads (default: 1)
-m,--mismatch <int>        number of mismatches to allow in each flanking region (defult: 0). An alignment is reported if there is a unique best alignment.
-s,--sam                   output aligned reads in .sam format
--rmdup                    remove PCR duplicates when reporting genotype counts
--no-anonymizer            don't run the anonymization step

Advanced options:
--min-read-length    minimum read length to process
--max-read-length    maximum read length to process
--error-rate         sequencing error rate to introduce
--fft-window-size    size of fft window (default: 24)
--fft-window-step    step size of sliding window (default: 12)
--lobe-threshold     threshold score to call a window periodic (defualt: 3)
--extend             length of flanking regions in the genome to align against (default: 100)
--minflank           minimum length of flanking region to try to align (default: 10)
--maxflank           length to trim flanking regions to if they exceed that length (default: 1000)
--max-diff-ref       maximum difference in length from the reference sequence to allow for alignment (default 50) (will take the absolute value)


This program takes in raw reads, detects and aligns reads containing microsatellites, and anonymizes reads that can be mapped uniquely to the genome.
"""

import os
import sys
import getopt
import random
debug = False
detPath = ""
alnPath = ""

####################################################
try:
    opts, args = getopt.getopt(sys.argv[1:], "f:t:g:o:hvqp:m:s:Nd:", ["files=","table=","genome=","out=","help","verbose","fastq","threads=","mismatch=","sam","fft-window-size=","fft-window-step=","lobe-threshold=","extend=","minperiod=","maxperiod=","minflank=","maxflank=","max-diff-ref=","rmdup","debug", "databse=","no-anonymizer","min-read-length=","max-read-length=","error-rate="])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)
args = [item[0] for item in opts]

if not(("-d" in args or "--database" in args or "-N" in args) and ("-f" in args or "--files" in args) and ("-t" in args or "--table" in args) and ("-g" in args or "--genome" in args) and ("-o" in args or "--out" in args)):
    usage()
    if not ("-h" in args or "--help" in args):
        if not ("-f" in args or "--files" in args):
            print "*****Please specifiy file(s)*****"
        if not ("-t" in args or "--table" in args):
            print "*****Please specify an STR table*****"
        if not ("-g" in args or "--genome" in args):
            print "*****Please specify a reference genome*****"
        if not ("-o" in args or "--out" in args):
            print "*****Please specify an output prefix name*****"
        if not ("-d" in args or "--database" in args or "-N" in args):
            print "*****Please specify an allele database or -N****"
    sys.exit(2)

# initialize variables

# mandatory arguments
infiles = []
tableFile = ""
genomeFile = ""
prefix = ""
databaseFile = ""

# options
verbose = False
fastq = False
sam = False
threads = 1
mismatch = 0
fftwindow = 24
fftstep = 12
lobethresh = 3
extend = 100
minperiod = 1
maxperiod = 8
minflank = 10
maxflank = 1000
maxdiffref = 50
rmdup = False
mask_n = False
anonymize = True
min_read_length = 50
max_read_length = 1000
error_rate = 0

# get values from arguments
for o,a in opts:
    if o == "--debug": debug = True
    if o == "--no-anonymizer": anonymize = False
    if o == "-d" or o == "--database":
        databaseFile = a
    if o == "-N": mask_n = True
    if o == "-f" or o == "--files":
        infiles = a.split(",")
    if o == "-t" or o == "--table":
        tableFile = a
    if o == "-g" or o == "--genome":
        genomeFile = a
    if o == "-o" or o == "--out": 
        prefix = a
    if o == "-v" or o == "--verbose":
        verbose = True
    if o == "-q" or o == "--fastq":
        fastq = True
    if o == "-p" or o == "--threads":
        try:
            threads = int(a)
        except:
            print "ERROR: --threads must be an integer"
            sys.exit(1)
        if threads < 1 or threads > 40:
            print "ERROR: Number of threads specified must be between 1 and 40"
            sys.exit(1)
    if o == "--minperiod":
        try:
            minperiod = int(a)
        except:
            print "ERROR: --minperiod must be an integer"
            sys.exit(1)
        if minperiod < 1:
            print "ERROR: --minperiod must be >= 1"
            sys.exit(1)
    if o == "--maxperiod":
        try:
            maxperiod = int(a)
        except:
            print "ERROR: --maxperiod must be an integer"
            sys.exit(1)
        if maxperiod < 1:
            print "ERROR: --maxperiod must be >=1"
            sys.exit(1)
    if o == "--extend":
        try:
            extend = int(a)
        except:
            print "ERROR: --extend must be an integer"
            sys.exit(1)
        if extend <0:
            print "ERROR: --extend cannot be negative"
            sys.exit(1)
        if extend > 500:
            print "ERROR: --extend must be < 500"
            sys.exit(1)
    if o == "--fft-window-size":
        try:
            fftwindow = int(a)
        except:
            print "ERROR: --fft-window-size must be an integer"
            sys.exit(1)
    if o == "--fft-window-step":
        try:
            fftstep = int(a)
        except:
            print "ERROR: --fft-window-step must be an integer"
            sys.exit(1)
    if o == "--lobe-threshold":
        try:
            lobethresh = float(a)
        except:
            print "ERROR: lobe threshold must be a float"
            sys.exit(1)
        if lobethresh <=0:
            print "ERROR: lobe threshold must be >=0"
            sys.exit(1)
    if o == "-m" or o == "--mismatch":
        try:
            mismatch = int(a)
        except: 
            print "ERROR: --mismatch must be an integer"
            sys.exit(1)
        if mismatch < 0:
            print "ERROR: cannot have a negative value for --mismatch"
            sys.exit(1)
    if o == "--minflank":
        try:
            minflank = int(a)
        except: 
            print "ERROR: --minflank must be an integer"
            sys.exit(1)
        if minflank < 0:
            print "ERROR: --minflank must be >= 0"
            sys.exit(1)
    if o == "--maxflank":
        try:
            maxflank = int(a)
        except:
            print "ERROR: --maxflank must be an integer"
            sys.exit(1)
        if maxflank < 0:
            print "ERROR: --maxflank must be >=0"
    if o == "--max-diff-ref":
        try:
            maxdiffref = int(a)
        except:
            print "ERROR: --max-diff-ref must be an integer"
            sys.exit(1)
    if o == "--rmdup":
        rmdup = True
    if o == "-s" or o == "--sam":
        sam = True
    if o == "--min-read-length":
        min_read_length = int(a)
    if o == "--max-read-length":
        max_read_length = int(a)
    if o == "--error-rate":
        error_rate = float(a)
        if error_rate > 1 or error_rate < 0:
            print "ERROR: --error-rate must be between 0 and 1"
            sys.exit(1)
    if o == "-h" or o == "--help": 
        usage()
        sys.exit(2)

if minperiod > maxperiod:
    print "ERROR: minperiod must be <= maxperiod"
    sys.exit(1)
if minflank > maxflank:
    print "ERROR: minflank must be <= maxflank"
    sys.exit(1)

###################################################

def command(cmd):
    """ run a command..."""
    if debug:
        print cmd
    else:
        os.system(cmd)

def main():
    ########### STEP 1/2/3 detection/alignment/genotyping ########
    if verbose:
        print "STEP 1/2/3: STR detection/alignment/genotyping..."    

    # run lobSTR
    files = ",".join(infiles)
    if fastq: fq = "-q"
    else: fq = ""
    lobSTR_cmd = "lobSTR -f %s -o %s -t %s -g %s --fft-window-size %s --minflank %s --maxflank %s --minperiod %s --maxperiod %s --fft-window-step %s --lobe-threshold %s -p %s %s -m %s --extend %s --max-diff-ref %s --min-read-length %s --max-read-length %s" %(files, prefix, tableFile, genomeFile, fftwindow, minflank, maxflank, minperiod, maxperiod, fftstep,lobethresh,threads, fq, mismatch, extend, maxdiffref, min_read_length, max_read_length)
    if sam:
        lobSTR_cmd += " --sam"
    if rmdup:
        lobSTR_cmd += " --rmdup"
    command(lobSTR_cmd)
    alignmentFile = "%s.aligned.tab"%prefix
    if verbose:
        print "STR detection/alignment/genotyping completed..."
        
    ########### STEP 4 anonymization #################
    if anonymize:
        if verbose: print "Performing anonymization..."
        strangers_cmd = "strangers -f %s -a %s -o %s %s --error-rate %s"%(files, alignmentFile, prefix, fq, error_rate)
        if mask_n: strangers_cmd += " --mask"
        else: strangers_cmd += " -d %s"% databaseFile
        command(strangers_cmd)
    if verbose: print "Done!"
main()

#!/usr/bin/env python
"""
python GenerateLobSTRRefBed.py <fastafile> <outfile>
"""

# Copyright (C) 2016 Melissa Gymrek <mgymrek@mit.edu> Thomas Willems <twillems@mit.edu>
#
# This file is part of lobSTR.
#
# lobSTR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# lobSTR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with lobSTR.  If not, see <http://www.gnu.org/licenses/>.

#
# A validation suite to check concordance with gold standard
# capillary electrophoresis calls.
#

import collections
import glob
import os
import pyfasta
import sys
import tempfile
from subprocess import Popen, PIPE, STDOUT

######################
# TRF params - at some point could be made options
TRF_MIN_SCORE = 5
TRF_MATCH_WT = 2
TRF_MISMATCH_PEN = 7
TRF_INDEL_PEN = 7
TRF_P_MATCH = 80
TRF_P_INDEL = 10
TRF_MAX_PERIOD = 500

# For fixing TRF output
MAX_PERIOD = 6
PATTERN_INDEX = 13
PERIOD_INDEX = 2
NREPEAT_INDEX = 3
CONS_SIZE_INDEX = 4
START_INDEX = 0
STOP_INDEX = 1
SCORE_INDEX = 7
NUM_TOKENS = 15

# For parsing TRF output
PARSE_NUM_TOKENS = 15
PARSE_SCORE_INDEX = 7
PARSE_PERIOD_INDEX = 4
PARSE_START_INDEX = 0
PARSE_STOP_INDEX = 1
PARSE_NREPEAT_INDEX = 3

######################

def PROGRESS(msg):
    sys.stderr.write(msg.strip() + "\n")

def ERROR(msg):
    sys.stderr.write(msg.strip() + "\n")
    sys.exit(1)

def RunCommand(cmd):
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, \
                  stderr=STDOUT, close_fds=True)
    ex = p.wait()
    if ex != 0:
        stdout, stderr = "", ""
        if p.stdout is not None: stdout = p.stdout.read()
        if p.stderr is not None: stderr = p.stderr.read()
    return ex

def CheckProgram(program):
    """ Check whether a program is installed """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        exe_file = os.path.join(path, program)
        if is_exe(exe_file): return True
    return False


# Create a new file containing the same information as contained in INPUT_FILE, 
# except that the period, pattern and number of repeats for entries 
# with incorrect patterns are corrected. Also filters out entries with a corrected
# period greater than MAX_PERIOD
# USE FOR FILES CONTAINING LONG PERIODS
def correct_pattern_errors(input_file, output_prefix, max_period=6):
    data    = open(input_file, "r")
 
    fail_count  = 0
    skip_count  = 0

    current_chrom = None
    output = None
    for line in data:
        if line.startswith("Sequence"):
            current_chrom = line.split()[1].strip()
            if output is not None: output.close()
            output = open(output_prefix + ".%s"%current_chrom, "w")
            continue

        tokens = line.strip().split()
        
        if len(tokens) != NUM_TOKENS:
            skip_count = skip_count + 1
            continue

        pattern = tokens[PATTERN_INDEX]
        subseq  = has_subseq_long(pattern)
        if subseq:
            new_patt = subseq
            tokens[NREPEAT_INDEX]   = str(float(tokens[NREPEAT_INDEX])*len(pattern)/len(new_patt))
            tokens[PATTERN_INDEX]   = new_patt
            tokens[PERIOD_INDEX]    = str(len(new_patt))
            tokens[CONS_SIZE_INDEX] = str(len(new_patt))
            new_line = ' '.join(tokens) + '\n'
            
            if len(new_patt) <= max_period:
                output.write(new_line)
            else:
                fail_count += 1
        else:
            if int(tokens[CONS_SIZE_INDEX]) <= max_period:
                output.write(line)
            else:
                fail_count += 1

    data.close()
    output.close()
    PROGRESS("# REMOVED PATTERNS = " + str(fail_count))

# Run the Z-algorithm on the provided sequence. Returns
# an array of length len(NMER), where arr[i] corresponds
# to the length of the sequence starting at position
# i that is a prefix of NMER
def zalgorithm(nmer):
    res = len(nmer)*[0]
    l   = -1
    r   = -1

    for i in xrange(1, len(nmer), 1):
        if i > r:
            j = i
            while j < len(nmer) and nmer[j] == nmer[j-i]:
                j = j + 1

            res[i] = j-i
            l      = i
            r      = j-1
        else:
            i_prime = i-l
            beta    = r-i+1

            if res[i_prime] < beta:
                res[i] = res[i_prime]
            elif res[i_prime] == beta:
                j = r + 1
                while j < len(nmer) and nmer[j] == nmer[j-i]:
                    j = j + 1

                res[i] = j-i
                l      = i
                r      = j-1                
            else:
                res[i] = beta
    res[0] = len(nmer)
    return res

# Returns a 2D array of length n, where
# the ith index contains all nmers of length i
def generate_nmers(max_n):
    nmers = []
    for i in xrange(max_n+1):
        nmers.append([])

    nmers[0].append('')
    for i in xrange(1, max_n+1, 1):
        for nmer in nmers[i-1]:
            nmers[i].append(nmer+'A')
            nmers[i].append(nmer+'C')
            nmers[i].append(nmer+'G')
            nmers[i].append(nmer+'T')
    return nmers

# Returns the substring from start to stop 
# (inclusive) of the provided sequence
def get_string(sequence, start, stop):
    if start > stop:
        return  sequence[start:]+sequence[:stop+1]
    else:
        return sequence[start:stop+1]


# Returns the repeating subunit of the provided sequence if it 
# consists of exactly 2 or more copies. Otherwise, returns False
# USE FOR LONG SEQUENCES
def has_subseq_long(nmer):
    prefix_lengths = zalgorithm(nmer)

    for k in xrange(1, len(nmer), 1):
        if len(nmer) % k == 0:
            match = True

            for segment in xrange(len(nmer)/k):
                coord = segment*k
               
                if prefix_lengths[coord] < k:
                    match = False
                    break

            if match:
                return nmer[0:k]
    return False


# Returns the repeating subunit of the provided sequence if it consits
# of exactly 2 or more exact copies. Otherwise, returns False
# ONLY USE FOR SHORT SEQUENCES 
def has_subseq(nmer):
    for k in xrange(1, len(nmer), 1):
        if len(nmer) % k == 0:
            seg_len = k
            
            seq = get_string(nmer, 0, (-1 + seg_len)%len(nmer))
            match = True

            for segment in xrange(len(nmer)/k):
                p1 = (segment*seg_len) % len(nmer)
                p2 = (-1 + (segment+1)*seg_len) % len(nmer)

                if get_string(nmer, p1, p2) != seq:
                    match = False
                    break

            if match:
                return seq
    return False


# Returns a dictionary that contains all kmers of length <= MAX_N
# that consist of a repeating subunit as keys and their smallest subunits as values 
def create_subseq_dict(max_n):
    nmers = generate_nmers(max_n)
    subseqs = {}
    for i in xrange(1, max_n+1, 1):
        for nmer in nmers[i]:
            subseq = has_subseq(nmer)
            if subseq:
                subseqs[nmer] = subseq
    return subseqs

# Create a new file containing the same information as contained in INPUT_FILE, 
# except that the period, pattern and number of repeats for entries 
# with incorrect patterns are corrected
def process_file(input_file, output_file, subseqs):
    data    = open(input_file, "r")
    output  = open(output_file, "w") 
    fail_count = 0

    for line in data:
        tokens = line.strip().split()
        pattern = tokens[PATTERN_INDEX]
        if pattern in subseqs:
            new_patt = subseqs[pattern]
            tokens[NREPEAT_INDEX]   = str(float(tokens[NREPEAT_INDEX])*len(pattern)/len(new_patt))
            tokens[PATTERN_INDEX]   = new_patt
            tokens[PERIOD_INDEX]    = str(len(new_patt))
            tokens[CONS_SIZE_INDEX] = str(len(new_patt))
            new_line = '\t'.join(tokens) + '\n'
            output.write(new_line)
            fail_count = fail_count + 1
        else:
            output.write(line)
    data.close()
    output.close()
    PROGRESS("# PATTERN CORRECTIONS = " + str(fail_count))

# Create a new file containing the same information as contained in INPUT_FILE, 
# except that some overlapping TR intervals are removed to get rid of all overlaps.
# Scanning from low to high coordinates on a chromosome, when an overlap occurs,
# the region with the higher score is kept while the region with the lower
# score is discarded and not written to the output file. If a discard file
# is provided, the removed regions are written to the discard file.
def resolve_overlaps(input_file, output_file, discard_file=None):
    data   = open(input_file, "r")
    lines  = []
    starts = []
    stops  = []
    scores = []
    for line in data:
        tokens = line.strip().split()
        starts.append(int(tokens[START_INDEX]))
        stops.append(int(tokens[STOP_INDEX]))
        scores.append(int(tokens[SCORE_INDEX]))
        lines.append(line)
    data.close()

    indices  = sorted(range(len(starts)), key = lambda x: (starts[x], stops[x]))
    max_stop = -1
    index    = -1
    score    = -1
    num_rem  = 0
    inc_tr   = len(starts)*[True]

    for i in xrange(len(starts)):
        if starts[indices[i]] <= max_stop:
            num_rem = num_rem + 1
            
            if scores[indices[i]] > score:
                max_stop      = stops[indices[i]]
                score         = scores[indices[i]]
                inc_tr[index] = False
                index         = indices[i]
            else:
                inc_tr[indices[i]] = False
        else:
            max_stop = stops[indices[i]]
            score    = scores[indices[i]]
            index    = indices[i]

    output = open(output_file, "w")
    if discard_file is not None:
        discard = open(discard_file, "w")

    for i in xrange(len(starts)):
        if inc_tr[i]:
            output.write(lines[i])
        elif discard_file is not None:
            discard.write(lines[i])

    output.close()
    if discard_file is not None:
        discard.close()

    PROGRESS("# REMOVED TANDEM REPEATS = " + str(num_rem))
    PROGRESS("TOTAL TANDEM REPEATS = " + str(len(starts)))

def parse_trf_file(input_file, store_vals=PARSE_NUM_TOKENS*[True]):
    data = open(input_file, "r")

    # start stop period num_repeats consensus_size percent_match percent_indel score a_frac c_frac g_frac t_frac entropy consensus sequence
    operations = [int, int, int, float, int, int, int, int, int, int, int, int, float, str, str]
    vals       = [[],  [],  [],  [],    [],  [],  [],  [],  [],  [],  [],  [],  [],    [],  []]
  
    skip_count = 0
    line_count = 0
    for line in data:
        tokens = line.strip().split()
        line_count = line_count + 1
        
        if len(tokens) != PARSE_NUM_TOKENS:
            skip_count = skip_count + 1
            continue

        for i in xrange(PARSE_NUM_TOKENS):
            if store_vals[i]:
                vals[i].append(operations[i](tokens[i])) 

    return_vals = []
    for i in xrange(PARSE_NUM_TOKENS):
        if store_vals[i]:
            return_vals.append(vals[i])

    PROGRESS("Skipped " + str(skip_count)+ "/" + str(line_count) + " lines")
    return return_vals


def create_filtered_trf_bed_file(input_prefix, output_file):
    max_period        = 6
    period_vals       = [1,  2,  3,  4,  5,  6]
    period_thresholds = [20, 22, 28, 28, 32, 34]

    # Get input files and chromosomes, store in input_files and chromosomes
    input_files = glob.glob(input_prefix + "*")
    chromosomes = [os.path.basename(f).split(".")[-1] for f in input_files]

    # Output fields and their associated indices in the input array
    #               start, end, period, refcopy, period, start, end, score, percA, percC, percG, percT, entropy, repseq
    data_indices = [0,     1,   4,      3,       4,      0,     1,   7,     8,     9,     10,    11,    12,      13]

    vals = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
    chrs = []
    for i in xrange(len(input_files)):
        PROGRESS("Processing " + input_files[i])
        file_vals = parse_trf_file(input_files[i])
        
        for j in xrange(len(file_vals[0])):
            if (file_vals[PARSE_PERIOD_INDEX][j] <= 6) and (file_vals[PARSE_SCORE_INDEX][j] > period_thresholds[file_vals[PARSE_PERIOD_INDEX][j]-1]):
                for k in xrange(len(vals)):
                    vals[k].append(file_vals[k][j])
                chrs.append(chromosomes[i])

    indices = range(len(chrs))
    indices = sorted(indices, key = lambda x: (chrs[x], vals[PARSE_START_INDEX][x], vals[PARSE_STOP_INDEX][x]))

    output = open(output_file, "w")
    for i in xrange(len(indices)):
        line = str(chrs[indices[i]])
        for j in xrange(len(data_indices)):
            line = line + "\t" + str(vals[data_indices[j]][indices[i]])
        line = line + "\n"
        output.write(line)
    output.close()

def analyze_overlaps(input_file, pass_file, fail_file):
    input_fh = open(input_file, "r")
    pass_fh = open(pass_file, "w")
    fail_fh = open(fail_file, "w")

    pass_count = 0
    fail_count = 0

    regions   = []
    min_start = 0
    max_stop  = -1
    for line in input_fh:
        tokens = line.strip().split("\t")
        chrom  = tokens[0]
        start, stop, period = map(int, tokens[1:4])
        num_repeats = float(tokens[4])
        score = float(tokens[8])
        if start <= max_stop:
            max_stop = max(max_stop, stop)
        else:
            # Analyze the previous set of overlapping regions
            region_size = max_stop - min_start + 1
            max_score   = 0
            max_index   = -1
            cov_frac    = 0.0
            for index,region in enumerate(regions):
                if region[4] >= max_score:
                    max_index = index
                    max_score = region[5]
                    cov_frac  = 1.0*(region[2]-region[1]+1)/region_size

            if max_index != -1:
                region = regions[max_index]
                if cov_frac > 0.85:
                    pass_fh.write("%s\t%s\t%d\t%d\t%.1f\n"%(region[0], region[1], region[2], region[3], region[4]))
                    pass_count += (len(regions) > 1)
                else:
                    for region in regions:
                        fail_fh.write("%s\t%s\t%d\t%d\t%.1f\n"%(region[0], region[1], region[2], region[3], region[4]))
                    fail_count += (len(regions) > 1)
            regions   = []
            min_start = start
            max_stop  = stop
        regions.append((chrom, start, stop, period, num_repeats, score))

    input_fh.close()
    pass_fh.close()
    fail_fh.close()
    PROGRESS("Pass: %s, Fail: %s"%(pass_count, fail_count))

def min_perm(seq):
    min_perm = seq
    for i in xrange(len(seq)):
        other = seq[i:]+seq[0:i]
        if other < min_perm:
            min_perm = other
    return min_perm

def longest_perfect_array(allele, motif):
    max_start  = -1
    max_length = -1
    for i in xrange(len(allele)):
        for offset in xrange(len(motif)):
            seq_index = i
            while seq_index < len(allele) and allele[seq_index] == motif[(offset+seq_index-i)%len(motif)]:
                seq_index += 1
            if seq_index - i > max_length:
                max_start  = i
                max_length = seq_index - i
    return max_length, max_start

# Complements for each of the 4 DNA bases
complement_table = {
    'A':'T',
    'C':'G',
    'G':'C',
    'T':'A',
    }

def complement(dna_sequence):
    res = ""
    for i in xrange(len(dna_sequence)):
        res = res + complement_table[dna_sequence[i]]
    return res

def reverse_complement(dna_sequence):
    return complement(dna_sequence[::-1])

def calc_alignment_purity(sequence, motif):
    """
    returns purity_score, num_errors, trf_score
    """
    # Traceback constants
    LEFT = -1
    DIAG = 0
    UP   = 1
    
    # Scoring constants
    MATCH    = 0
    MISMATCH = 1
    INDEL    = 1

    # Create a "perfect" sequence consisting of copies of the motif
    perfect_seq  = (int)(len(sequence)*3.0/len(motif) + 1)*motif
        
    # Init trf score
    trf_score = 0

    # Create scoring and traceback matrices
    score_matrix = []
    trace_matrix = []
    for i in xrange(len(sequence)+1):
        score_matrix.append((len(perfect_seq)+1)*[None])
        trace_matrix.append((len(perfect_seq)+1)*[None])
    
    # Initialize matrices with no endspace penalty for the "perfect" sequence
    score_matrix[0][0] = 0
    trace_matrix[0][0] = None
    for j in xrange(1, len(perfect_seq)+1):
        score_matrix[0][j] = 0
        trace_matrix[0][j] = LEFT
    for i in xrange(1, len(sequence)+1):
        score_matrix[i][0] = i*INDEL
        trace_matrix[i][0] = UP

    # Perform alignment
    for i in xrange(1, len(sequence)+1):
        for j in xrange(1, len(perfect_seq)+1):
            v1 = score_matrix[i-1][j]+INDEL
            v2 = score_matrix[i][j-1]+INDEL
            v3 = score_matrix[i-1][j-1] + (MATCH if sequence[i-1] == perfect_seq[j-1] else MISMATCH)

            if v1 < v2:
                if v1 < v3:
                    score_matrix[i][j] = v1
                    trace_matrix[i][j] = UP
                else:
                    score_matrix[i][j] = v3
                    trace_matrix[i][j] = DIAG
            else:
                if v2 < v3:
                    score_matrix[i][j] = v2
                    trace_matrix[i][j] = LEFT
                else:
                    score_matrix[i][j] = v3
                    trace_matrix[i][j] = DIAG
    
    # Find the optimal alignment that uses the entire provided sequence
    min_score = score_matrix[len(sequence)][0]
    min_index = 0
    for j in xrange(len(perfect_seq)+1):
        if score_matrix[len(sequence)][j] < min_score:
            min_score = score_matrix[len(sequence)][j]
            min_index = j          
            
    # Traceback the optimal alignment
    pseq_alignment = ""
    qseq_alignment = ""
    i = len(sequence)
    j = min_index   
    while i > 0:
        if trace_matrix[i][j] == LEFT:
            pseq_alignment += perfect_seq[j-1]
            qseq_alignment += "-"
            j -= 1
            trf_score = trf_score - TRF_INDEL_PEN
        elif trace_matrix[i][j] == UP:
            pseq_alignment += "-"
            qseq_alignment += sequence[i-1]
            i -= 1
            trf_score = trf_score - TRF_INDEL_PEN
        elif trace_matrix[i][j] == DIAG:
            pseq_alignment += perfect_seq[j-1]
            qseq_alignment += sequence[i-1]
            i -= 1
            j -= 1
            if pseq_alignment[-1] == qseq_alignment[-1]:
                trf_score = trf_score + TRF_MATCH_WT
            else: trf_score = trf_score - TRF_MISMATCH_PEN
        else:
            ERROR("Invalid condition encountered in alignment traceback. Exiting...")
    pseq_alignment = pseq_alignment[::-1]
    qseq_alignment = qseq_alignment[::-1]
    num_errors     = min_score
    purity_score   = 1.0 - 1.0*num_errors/(len(sequence))
    return purity_score, num_errors, trf_score

def annotate_reference(reffasta, input_file, output_file):
    input_fh = open(input_file, "r")
    output_fh = open(output_file, "w")
    fasta = pyfasta.Fasta(reffasta)
    keymap = {}
    for k in fasta.keys():
        keymap[k.split()[0]] = k

    for line in input_fh:
        tokens      = line.strip().split()
        chrom       = tokens[0]
        start, stop = map(int, tokens[1:3])
        period      = int(tokens[3])
        seq = fasta[keymap[chrom]][start-1:stop].upper()        
        motif_counts = collections.defaultdict(int)
        for i in xrange(len(seq)-period):
            motif_counts[min_perm(seq[i:i+period])] += 1
        most_common_motif    = sorted(motif_counts.items(), key = lambda x: x[1])[-1][0]
        trf_score            = calc_alignment_purity(seq, most_common_motif)[2]
        pure_len             = longest_perfect_array(seq, most_common_motif)[0]
        most_common_motif    = min(most_common_motif, min_perm(reverse_complement(most_common_motif)))
        output_fh.write(line.strip() + "\t" + str(most_common_motif) + "\t" + str(trf_score) + "\t" + str(pure_len) + "\n")
    input_fh.close()
    output_fh.close()

############################
def main():
    try:
        FASTAFILE = sys.argv[1]
        OUTFILE = os.path.abspath(sys.argv[2])
    except:
        print __doc__
        sys.exit(1)

    # Check for required programs
    if not CheckProgram("trf"):
        ERROR("Could not find trf installed on the PATH")
    # Check input file exists
    if not os.path.exists(FASTAFILE):
        ERROR("%s does not exist"%FASTAFILE)
    name = os.path.basename(FASTAFILE)

    # Set up temporary directory
    tmpdir = tempfile.mkdtemp(prefix="lobref.")
    os.chdir(tmpdir)

    # Run TRF
    PROGRESS("Running TRF on %s"%FASTAFILE)
    os.mkdir("source_results")
    os.chdir("source_results")
    trf_cmd = "trf %s %s %s %s %s %s %s %s -h -d"%(FASTAFILE, TRF_MATCH_WT, TRF_MISMATCH_PEN, TRF_INDEL_PEN, \
                                                 TRF_P_MATCH, TRF_P_INDEL, TRF_MIN_SCORE, TRF_MAX_PERIOD)
    if RunCommand(trf_cmd):
        #ERROR("Error running %s"%trf_cmd)
        PROGRESS("WARNING: trf exited with code 1")
    mv_trf_cmd = "mv %s.%s.%s.%s.%s.%s.%s.%s.dat %s"%(name, TRF_MATCH_WT, TRF_MISMATCH_PEN, TRF_INDEL_PEN, \
                                                      TRF_P_MATCH, TRF_P_INDEL, TRF_MIN_SCORE, TRF_MAX_PERIOD, name)
    if RunCommand(mv_trf_cmd):
        ERROR("Error running %s"%mv_trf_cmd)
    os.chdir("..")

    # Fix errors in TRF patterns and remove loci with P > 6
    PROGRESS("Adjusting TRF patterns")
    os.mkdir("fixed_results")
    correct_pattern_errors(os.path.join(tmpdir, "source_results", name), \
                               os.path.join(tmpdir, "fixed_results", name), \
                               max_period=MAX_PERIOD)

    # Filter regions by scoring threshold
    PROGRESS("Filter regions by scoring threshold")
    os.mkdir("filtered_results")
    create_filtered_trf_bed_file(os.path.join(tmpdir, "fixed_results", name), \
                                     os.path.join(tmpdir, "filtered_results", "all_chrs"))


    # Split back up by chromosome
    PROGRESS("Split back up by chromosome")
    split_cmd = "awk '{print > \"%s/\"$1}' %s/all_chrs"%(os.path.join(tmpdir, "filtered_results"), \
                                                             os.path.join(tmpdir, "filtered_results"))
    if RunCommand(split_cmd):
        ERROR("Error running %s"%split_cmd)

    # Combine overlapping TRFs contingent upon mering criteria
    PROGRESS("Combine overlapping TRFs")
    os.mkdir("merged_results")
    for chrfile in glob.glob(os.path.join(tmpdir, "filtered_results/*")):
        if "all_chrs" in chrfile: continue
        analyze_overlaps(chrfile, \
                             os.path.join(tmpdir, "merged_results", os.path.basename(chrfile)+".pass"), \
                             os.path.join(tmpdir, "merged_results", os.path.basename(chrfile)+".fail"))

    # Create final reference
    PROGRESS("Create final reference")
    merge_cmd = "cat %s/*.pass | sort -k 1,1V -k 2,2n > %s"%(os.path.join(tmpdir, "merged_results"), \
                                                                 os.path.join(tmpdir, "passing_str_reference.bed"))
    if RunCommand(merge_cmd):
        ERROR("Error running %s"%merge_cmd)

    # Annotate with motif and number of imperfections
    annotate_reference(FASTAFILE, os.path.join(tmpdir, "passing_str_reference.bed"), os.path.join(tmpdir, "ann_ref.bed"))

    # Generate lobSTR table
    # TODO make score
    col_cmd = "cat %s | awk '{print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t\" $4 \"\t\" $2 \"\t\" $3 \"\t\" $7 \"\tNA\tNA\tNA\tNA\tNA\t\" $6 \"\t.\"}' > %s"%(os.path.join(tmpdir, "ann_ref.bed"), OUTFILE)
    if RunCommand(col_cmd):
        ERROR("Error running %s"%col_cmd)

if __name__ == "__main__":
    main()

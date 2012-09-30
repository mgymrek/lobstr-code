def usage():
    print """
python lobSTR_allelotype_checks.py -f <genotypes.tab file>  [--plot] [--mincov <INT>]

Note: this script only works on files generated with lobSTR v2.0.0 or later.

Runs the following lobSTR quality checks on output of lobSTR allelotyping:
0. Number of genotype calls (cov >= 1)
1. Number of genotype calls (cov >= 5)
2. Mean coverage (only fully spanning reads)
3. Mean % agreeing reads
4. STR period vs. fraction in each category (0/0, 0/1, 1/1, 1/2)

If the --plot option is specified, the following plots are generated:
0. Coverage plot
1. STR period vs. fraction in each category (0/0, 0/1, 1/1, 1/2)
2. Distribution of distance between alleles (by period)
3. Distribution of percentage of agreeing reads
4. Distribution of allelotype scores
5. Distribution of diff from ref of alleles
6. Non-unit repeat alleles plot
7. Stutter noise plots

Prints all results to stdout
Note stutter is determined using high confidence homozygous calls

"""

###########################

from checks_utils import *
import math

###########################
MIN_COV = 5
MOTIF_COL = 3
PERIOD_COL = 5
REF_COL = 6
ALLELES_COL = 7
COV_COL = 8
AGREE_COL = 9
ALL_READS_COL = 11
###########################

try:
    opts, args = getopt.getopt(sys.argv[1:], "hvf:", ["help","verbose","mincov=","plot","debug"])
except getopt.GetoptError, err:
    print str(err)
    usage()
    sys.exit(2)

args = [item[0] for item in opts]

if ((not "-f" in args)):
    usage()
    sys.exit(2)

# initialize variables
verbose = False
filename = ""
plot = False
debug = False

# set variables
for o,a in opts:
    if o == "-f": filename = a
    if o == "--mincov": MIN_COV = int(a)
    if o == "--plot": plot = True
    if o == "--debug": debug = True
    if o == "-v" or o == "--verbose": verbose = True
    if o == "--help" or o == "-h":
        usage()
        sys.exit(0)

###########################

if plot:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams['pdf.fonttype'] = 42
###########################

def GetMeanAgreeingReads(filename):
    """
    Get mean % of agreeing reads
    """
    cmd = "grep -v version %s | awk \'($%s > 0)\' | "\
        "awk \'{print $%s/$%s}\' | awk \'{mean += $1}END{print mean/NR}\'"\
        %(filename, COV_COL, AGREE_COL, COV_COL)
    return ExecuteCmd(cmd, debug)

def GetNumGenotypesByCov(filename, cov):
    """
    Get number of STRs covered by $cov or more reads
    """
    cmd = "cat %s | grep -v version | awk \'($%s >=%s)\' | wc -l"%(filename, COV_COL, cov)
    return ExecuteCmd(cmd, debug)


def ClassifyAllelotypes(filename):
    """
    Get allelotype categories for each period
    0: homozygous ref
    1: heterozygous nref/ref
    2: homozygous nref/nref
    3: heterozygous nref/nref
    """
    sample = filename
    f = open(filename, "r")
    line = f.readline()
    while "version" in line: line = f.readline()
    # period -> type -> count
    gendict = {2:[0,0,0,0], 3:[0,0,0,0], 4:[0,0,0,0],\
                   5:[0,0,0,0], 6:[0,0,0,0]}
    # period -> modunit -> count
    moddict = {2:[0,0,0,0,0,0], 3:[0,0,0,0,0,0], 4:[0,0,0,0,0,0],\
                5:[0,0,0,0,0,0], 6:[0,0,0,0,0,0]}
    # list of [period, motif, stutter(0/1), steps, length]
    stutter_info_list = []
    while line != "":
        items = line.strip().split("\t")
        motif = items[MOTIF_COL-1]
        period = int(items[PERIOD_COL-1])
        coverage = int(items[COV_COL-1])
        agree = int(items[AGREE_COL-1])
        if coverage > 0:
            percagree = agree*1.0/coverage
        else:
            percagree = 0
        ref = float(items[REF_COL-1])
        try:
            if coverage >= MIN_COV:
                alleles = [int(item) for item in items[ALLELES_COL-1].split(",")]
            else: alleles = []
        except:
            # partially spanning - don't process
            alleles = []
        # get mod info
        for item in alleles:
            moddict[period][item%period] = moddict[period][item%period]+1
        # get stutter info (confident homozygous calls)
        if len(alleles) == 2:
            if ( coverage >= MIN_COV and alleles[0] == alleles[1] and percagree >= 0.9):
                all_reads = items[ALL_READS_COL-1].split("/")
                genotype = int(items[ALLELES_COL-1].split(",")[0])
                for allele in all_reads:
                    call = int(allele.split(":")[0])
                    numreads = int(allele.split(":")[1])
                    step = (call-genotype)*1.0/period
                    length = ref*period + call
                    if call == genotype:
                        stutter = 0
                    else: stutter = 1
                    if abs(step)*1.0 <=3:
                        for i in range(numreads):
                            stutter_info_list.append([period, motif, stutter, step, int(length)])
        if alleles != []:
            # homozygous
            if alleles[0] == alleles[1]:
                # Case 0: homozygous ref
                if alleles[0] == 0:
                    gendict[period][0] = gendict[period][0] + 1
                # Case 2: homozygous nref
                else:
                    gendict[period][2] = gendict[period][2] + 1
            # heterozygous
            else:
                # Case 1: heterozygous nref/ref
                if alleles[0] == 0 or alleles[1] == 0:
                    gendict[period][1] = gendict[period][1] + 1
                # Case 3: heterozygous nref/nref
                else:
                    gendict[period][3] = gendict[period][3] + 1
        line = f.readline()
    f.close()
    # generate results string
    results_string = "%s\tPeriod\t0/0\t0/1\t1/1\t1/2\n"%sample
    for period in range(2,7):
        loci = sum([gendict[period][i] for i in range(4)])
        try:
            percs = [gendict[period][i]*1.0/loci for i in range(4)]
        except:
            percs = [-1 for i in range(4)]
        results_string += "%s\tPeriod:%s\t%s\t%s\t%s\t%s\n"\
            %(sample,period, percs[0], percs[1], percs[2], percs[3])
    return results_string, gendict, moddict, stutter_info_list

def GetMeanCoverage(filename, zipped=False):
    """
    Get mean coverage of an STR (non-partial only)
    """
    if zipped:
        cmd = "zcat %s | grep -v version | awk '($8 > 0)'|"\
            "  awk \'{matesum += $(8)}END"\
            "{print matesum/NR}\'"%filename 
    else:
        cmd = "cat %s | grep -v version | awk '($8 > 0)' |"\
            "  awk \'{matesum += $(8)}END"\
            "{print matesum/NR}\'"%filename 
    return ExecuteCmd(cmd, debug)

def GetCoverages(filename):
    """
    Get coverage values for STRs
    """
    cmd = "cat %s | grep -v version | awk '($8 != 0)' | cut -f 8"%filename
    cov_stdout = ExecuteCmd(cmd, debug)
    covs = [int(item) for item in GetListFromStdout(cov_stdout)]
    return covs

def GetAlleleDiffs(filename, period):
    """
    Get distance between alleles (in bp) at het. loci (by period)
    """
    cmd = "cat %s | grep -v version | awk '($8 >= %s)' | awk '($5 == %s)' | cut -f 7 | "\
        "sed 's/,/\\t/g' | awk '($1!=$2)' | awk '(x=$2-$1) {print (x<0?-1*x:x)}'"%(filename, MIN_COV, period)
    diffs_stdout = ExecuteCmd(cmd, debug)
    try:
        diffs = [int(item) for item in GetListFromStdout(diffs_stdout)]
    except: diffs = []
    return diffs

def GetPercAgreeingList(filename):
    cmd = "grep -v version %s | awk \'($8 > 0)\' | "\
        "awk \'{print $9/$8}\' "\
        %(filename)
    agreeing_stdout = ExecuteCmd(cmd, debug)
    agreeing = [float(item) for item in GetListFromStdout(agreeing_stdout)]
    return agreeing

def GetAllelotypingScores(filename):
    cmd = "grep -v version %s | awk \'($8 > 0)\' | "\
        "cut -f 12 | grep -v inf "\
        %(filename)
    scores_stdout = ExecuteCmd(cmd, debug)
    scores = [float(item) for item in GetListFromStdout(scores_stdout)]
    return scores

def GetDiffFromRef(filename):
    cmd = "grep -v version %s | awk \'($8 > 0)\' | "\
        "cut -f 7 | grep -v chrX | grep -v chrY | grep -v NA | sed 's/,/\\n/g' "\
        %(filename)
    diffs_stdout = ExecuteCmd(cmd, debug)
    diffs = [float(item) for item in GetListFromStdout(diffs_stdout)]
    return diffs

###########################
def main():
    lobstr_genotype_file = filename

    # 0. Number of genotype calls (cov >= 1)
    num_total = GetLineCount(lobstr_genotype_file, debug)
    if not debug:
        print("%s\tNumber STR covered at least once (includes partially spanning)\t%s"%(lobstr_genotype_file, num_total))

    # 1. Number of genotype calls (cov >= 1, 5, 10)
    cov_levels = [1, 5, 10]
    for cov in cov_levels:
        num_genotype_calls = GetNumGenotypesByCov(lobstr_genotype_file, cov)
        if not debug:
            print("%s\tNumber STR loci coverage >= %s\t%s"%(lobstr_genotype_file, cov, num_genotype_calls))

    # 2. Mean coverage
    mean_cov = GetMeanCoverage(lobstr_genotype_file)
    if not debug:
        print("%s\tMean STR coverage (only non-partial reads counted)\t%s"%(lobstr_genotype_file, mean_cov))

    # 3. Mean % agreeing reads
    mean_agreeing = GetMeanAgreeingReads(lobstr_genotype_file)
    if not debug:
        print("%s\tMean percent agreeing reads\t%s"%(lobstr_genotype_file, mean_agreeing))

    # 4. STR period vs. fraction in each category (0/0, 0/1, 1/1, 1/2)
    print("%s\tSTR period vs. fraction in each category"%lobstr_genotype_file)
    if not debug:
        allelotype_classes, allelotypes_class_dict, moddict, stutter_info_list = ClassifyAllelotypes(lobstr_genotype_file)
        print(allelotype_classes.strip())

    if plot:
        # 0. Coverage plot
        plt.clf()
        covs = GetCoverages(lobstr_genotype_file)
        n, bins, patches = plt.hist(covs,  range=[0,50])
        plt.xlabel("Coverage")
        plt.ylabel("# STRs")
        plt.title("")
        plt.savefig("%s.coverage.pdf"%lobstr_genotype_file)

        # 1. STR period vs. fraction in each category (0/0, 0/1, 1/1, 1/2)
        plt.clf()
        types = range(4)
        periods = range(2,7)
        bottoms = [0]*5
        colors = {0:"red", 1:"blue", 2: "green", 3: "orange"}
        for i in types:
            type_counts = []
            for period in periods:
                counts = (allelotypes_class_dict.get(period,[0,0,0,0]))
                num = counts[i]
                if sum(counts) != 0:
                    num = num*1.0/sum(counts)
                type_counts.append(num)
            plt.bar(periods, type_counts, color = colors[i], bottom=bottoms,  align='center')
            bottoms = [type_counts[n] + bottoms[n] for n in range(len(type_counts))]
        plt.ylim(0,1)
        plt.xlabel("Period")
        plt.ylabel("Percentage of loci")
        plt.savefig("%s.locitypes.pdf"%lobstr_genotype_file)

        # 2. Distribution of distance between alleles (by period)
        colors = {2:"green", 3:"orange", 4:"red", 5:"blue", 6:"purple"}
        for period in xrange(2,7):
            plt.clf()
            allele_diffs = GetAlleleDiffs(lobstr_genotype_file, period)
            if (len(allele_diffs) > 0):
                n, bins, patches = plt.hist(allele_diffs, facecolor=colors[period])
                plt.xlabel("Distance between alleles (bp) (Period %s)"%period)
                plt.ylabel("# STRs")
                plt.title("")
                plt.savefig("%s.alleledist.per%s.pdf"%(lobstr_genotype_file, period))

        # 3. Distribution of perc agreeing
        plt.clf()
        perc_agreeing = [item for item in GetPercAgreeingList(lobstr_genotype_file) if item > 0]
        n, bins, patches = plt.hist(perc_agreeing, align='center')
        plt.xlabel("Percentage of reads agreeing")
        plt.ylabel("# STRs")
        plt.title("")
        plt.savefig("%s.percagree.pdf"%lobstr_genotype_file)

        # 4. Distribution of allelotype scores
        plt.clf()
        scores = [item for item in GetAllelotypingScores(lobstr_genotype_file) if item > 0]
        n, bins, patches = plt.hist(scores)
        plt.xlim(0,1)
        plt.xlabel("Allelotype scores")
        plt.ylabel("# STRs")
        plt.title("")
        plt.savefig("%s.allelotypescores.pdf"%lobstr_genotype_file)
        
        # 5. Distribution of diff from ref of alleles
        plt.clf()
        diffs = GetDiffFromRef(lobstr_genotype_file)
        n, bins, patches = plt.hist(diffs, range=[-50,50])
        plt.xlim(-50,50)
        plt.xlabel("Diff from ref")
        plt.ylabel("# STRs")
        plt.title("")
        plt.savefig("%s.difffromref.pdf"%lobstr_genotype_file)

        # 6. Non-unit repeat alleles plot
        plt.clf()
        colors = {2:"green", 3:"orange", 4:"red", 5:"blue", 6:"purple"}
        for p in xrange(2,7):
            y = [math.log10(item) for item in moddict[p] if item > 0]
            x = range(len(y))
            plt.plot(x, y, 'o', color=colors[p])
            plt.plot(x, y, color=colors[p], lw=2)
        plt.xlabel("Allele bp modulo motif period")
        plt.ylabel("log10 Number of loci")
        plt.title("")
        plt.savefig("%s.unit.pdf"%lobstr_genotype_file)

        # 7. Stutter noise
        # 7a. Stutter by period
        stutter_probs_by_period = []
        for period in xrange(2,7):
            all_reads = ([item for item in stutter_info_list if item[0] == period])
            stutter_reads = ([item for item in all_reads if item[2]])
            try:
                sprob = len(stutter_reads)*1.0/len(all_reads)
            except: sprob = 0
            stutter_probs_by_period.append(sprob)
        plt.clf()
        plt.bar(range(2,7), stutter_probs_by_period, align='center')
        plt.xlabel("Period")
        plt.ylabel("Stutter probability")
        plt.savefig("%s.stutter_by_period.pdf"%lobstr_genotype_file)

        stutter_probs_by_period = []
        for period in xrange(2,7):
            all_reads = ([item for item in stutter_info_list if item[0] == period and item[3] == int(item[3])])
            stutter_reads = ([item for item in all_reads if item[2]])
            try:
                sprob = len(stutter_reads)*1.0/len(all_reads)
            except: sprob = 0
            stutter_probs_by_period.append(sprob)
        plt.clf()
        plt.bar(range(2,7), stutter_probs_by_period, align='center')
        plt.xlabel("Period")
        plt.ylabel("Stutter probability")
        plt.savefig("%s.stutter_by_period.unit.pdf"%lobstr_genotype_file)

        # 7b. Stutter by length
        binsize = 5.0
        length_bins = range(25,80,int(binsize))
        stutter_probs_by_length = []
        for i in range(len(length_bins)):
            ldata = [item for item in stutter_info_list if item[4] > length_bins[i]-binsize/2 and item[4] <= length_bins[i]+binsize/2]
            stutter_data = [item for item in ldata if item[2]]
            try:
                sprob = len(stutter_data)*1.0/len(ldata)
            except: sprob = 0
            stutter_probs_by_length.append(sprob)
        plt.clf()
        plt.xlabel("Array length")
        plt.ylabel("Stutter probability")
        plt.plot(length_bins, stutter_probs_by_length, 'o')
        plt.savefig("%s.stutter_by_length.pdf"%lobstr_genotype_file)

        # 7c. Length of stutter
        all_stutter = [item for item in stutter_info_list if item[2]]
        stutter_steps = [round(item[3]) for item in all_stutter]
        plt.clf()
        plt.hist(stutter_steps, bins=range(-4,5), align="center")
        plt.xlabel("Stutter step sizes")
        plt.ylabel("Frequency")
        plt.title("")
        if (len(all_stutter) >0):
            plt.savefig("%s.stutter_lengths.pdf"%lobstr_genotype_file)
        
        # 7d. Heatmap of stutter by length for each period
        for period in xrange(2,7):
            # generate matrix of (true allele, allele seen)
            pdata = [item for item in stutter_info_list if item[0] == period]
            # get dimensions
            true_alleles = [item[4] for item in pdata]
            try:
                min_true = min(true_alleles)
                max_true = max(true_alleles)
            except:
                min_true = 0
                max_true = 0
            obs_alleles = [item[4]+item[3]*period for item in pdata]
            try:
                min_obs = min(obs_alleles)
                max_obs = max(obs_alleles)
                min_all = min([min_obs, min_true])
                max_all = max([max_obs, max_true])
            except:
                min_obs = 0
                max_obs = 0
                min_all = 0
                max_all = 0
            num_all = max_all-min_all+1
            data = [0 for i in xrange(int(num_all*num_all))]
            data = np.array(data)
            data.shape = (num_all, num_all)
            # fill in data
            for read in pdata:
                true = read[4]
                obs = read[4]+read[3]*period
                data[true-min_all,obs-min_all] = data[true-min_all,obs-min_all] + 1
            data = np.log(data+0.001)
            plt.clf()
            fig = plt.figure()
            plt.pcolormesh(data, cmap="Blues")
            plt.colorbar()
            if (len(pdata)>0):
                plt.savefig("%s.stutterheat.period%s.pdf"%(filename, period))

main()

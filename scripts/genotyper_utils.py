import math
import statsmodels.api as sm
from numpy import array
import pysam
import sys

# used to prevent problems coverting to log scale 
SMALL_CONST = 1e-100
# used as a guess for the maximal length. not critical
MAX_STR_LEN = 500
# other consts
CHR_Y = "chrY"
CHR_X = "chrX"
TRAIN_MAX_LEN = 25 
MINIMAL_OBS_FOR_TRUNC = 10
MIN_HET_FREQ = 0.20

######## methods ############
def GetTagValue(aligned_read, tag):
    for item in aligned_read.tags:
        if item[0] == tag:
            return item[1]
    if (tag == "XD"):
        print "No XD tag found. This bam file did not come from lobSTR. You must run convertBam.py first to add the necessary samtools flags"
        sys.exit(2)
    return None

def invLogit(x):
    """
    invLogit(x) returns the inverse of the logistic function.
    the logistic function logit(p) is given by log(p/(1-p))
    so the inverse is exp(x)/(1+exp(x))
    input: x, a number in R
    output: a probability
    """
    exp = math.exp(x)
    return( exp/ (1+exp))

def GetCopyNumbers(aligned_read_list):
    copy_numbers = []
    for aligned_read in aligned_read_list:
        try:
            copy_numbers.append(GetTagValue(aligned_read, "XC")+
                                GetTagValue(aligned_read, "XD")*1.0/
                                len(GetTagValue(aligned_read, "XR")))
        except: pass
    return copy_numbers

def ppois(x, mean):
    """
    ppois(x, mean) returns the probability that a Poisson distribution
    with the given mean would generate the value x
    input: x (the observed value)
           mean (the expected value, usually called lambda)
    output: A probability
    """
    assert(mean>0)
    assert(int(x)==x)
    if x<0: return 0
    p = math.exp(-mean)
    # implemented with a loop to avoid overflows
    for i in xrange(x):
        p *= mean
        p /= i+1
    return p


class ReadContainer:
    """ Keep track of aligned reads """
    def __init__(self):
        """ initializer does nothing """
        # map of (chrom, start) -> list of AlignedReads
        self.aligned_str_map_ = {}

    def AddReadsFromFile(self, bamfile):
        """ add each read from a bam file """
        samfile = pysam.Samfile( bamfile, "rb" )
        for read in samfile.fetch():
            chrom = samfile.getrname(read.rname)
            start = read.opt("XS")
            end = read.opt("XE")
            if int(end) > int(start):
                self.AddRead(chrom, start, end, read)

    def AddRead(self, chrom, start, end, aligned_read):
        """ read a line from the bam file and add to the genotyper """
        key = (chrom, start, end)
        if key in self.aligned_str_map_:
            self.aligned_str_map_[key].append(aligned_read)
        else:
            self.aligned_str_map_[key] = [aligned_read]

    def RemovePCRDuplicates(self):
        """ remove pcr duplicates from self.aligned_str_map_"""
        for str_locus in self.aligned_str_map_:
            pcr_duplicates = {}
            # group into duplicates
            for aligned_read in self.aligned_str_map_[str_locus]:
                key = (aligned_read.pos, aligned_read.qend)
                pcr_duplicates[key] = pcr_duplicates.get(key,[])+[aligned_read]
            # choose one rep from each group
            reads_after_rmdup = []
            for duplicate_group in pcr_duplicates:
                reads_after_rmdup.append(self.GetRepRead(pcr_duplicates[duplicate_group]))
            # reset that entry in the dictionary
            self.aligned_str_map_[str_locus] = reads_after_rmdup

    def GetRepRead(self,aligned_read_list):
        """
        return read with  majority copy num.
        Use quality scores to break ties
        """
        copy_number_to_reads = {}
        for aligned_read in aligned_read_list:
            diff = GetTagValue(aligned_read, "XD")
            copy_number_to_reads[diff] = copy_number_to_reads.get(diff,[])+[aligned_read]
        # majority ?
        copy_number_to_counts = [(len(copy_number_to_reads[copy_num]), copy_num) for copy_num in copy_number_to_reads]
        copy_number_to_counts.sort(reverse=True)
        top_count = copy_number_to_counts[0][0]
        top_copy_number = copy_number_to_counts[0][1]
        candidate_copy_numbers = [item for item in copy_number_to_counts if item[0] == top_count]
        if len(candidate_copy_numbers) == 1:
            return copy_number_to_reads[top_copy_number][0]
        # else use quality scores
        top_quality_score = 0
        top_copy_number = None
        for (count, diff) in candidate_copy_numbers:
            average_quality = self.GetAverageQuality(copy_number_to_reads[diff])
            if average_quality > top_quality_score:
                top_quality_score = average_quality
                top_copy_number = diff
        return copy_number_to_reads[top_copy_number][0]

    def GetAverageQuality(self, read_list):
        """ return average quality score of all reads in the list """
        if len(read_list) == 0: return 0
        total_quality = 0
        for read in read_list:
            qual_score = self.GetScore(read.qual)
            total_quality += qual_score
        return total_quality*1.0/len(read_list)
            
    def GetScore(self, quality_string):
        """ get quality score for a single read """
        if len(quality_string) == 0: return 0
        total_quality = 0
        for pos in quality_string:
            score = ord(pos) - 33
            total_quality += score
        return total_quality*1.0/len(quality_string)

class NoiseModel:
    """
    noiseModel: fit a noise model from a list of files or read the parameters
    from a file. 
    """

    def __init__(self, read_container = None):
        self.read_container_ = ReadContainer()
        if (read_container != None):
            self.read_container_ = read_container

    def train(self):
        """
        estimate the parameters of the noise model
        """
        self.chrs = []
        self.unitSize = []
        self.reads = []

        # populate chrs, unitSize, and reads from the read_container
        for str_region in self.read_container_.aligned_str_map_.keys():
            self.chrs.append(str_region[0])
            # get unitsize from tag XR length in first entry
            self.unitSize.append(len(self.read_container_.aligned_str_map_[str_region][0].opt("XR")))
            # get reads from diff from ref for each
            copy_numbers = GetCopyNumbers(self.read_container_.aligned_str_map_[str_region])
            for aligned_read in self.read_container_.aligned_str_map_[str_region]:
                copy_numbers.append(aligned_read.opt("XC")+
                                    aligned_read.opt("XD")*1.0/
                                    len(aligned_read.opt("XR")))
            self.reads.append(copy_numbers)
        # get the reads from chrX and chrY with a unique mode
        XY_indices = [x for x in range(len(self.chrs)) if self.chrs[x] in [CHR_Y,CHR_X]\
                      and self.hasUniqueMode(self.reads[x])!= False]
        chrD = {}
        for i in XY_indices: chrD[self.chrs[i]] = chrD.get(self.chrs[i],0)+1
        
        # create data for training
        # training data format:
        # [mode, period, copy number]
        trainingData = []
        modeD = {}
        for XY_index in XY_indices:
            mode = math.floor(self.hasUniqueMode(self.reads[XY_index]))
            for read in self.reads[XY_index]:
                trainingData.append([mode,self.unitSize[XY_index],math.floor(read)])
        self.trainingData = trainingData
        self.fitMutProb() # probability to mutate based on STR unit length
        #self.fitDecProb() # think removing this step
        self.fitPois()
        
    def hasUniqueMode(self,readSet):
        """
        gets a list of reads and returns the mode if it is unique and F otherwise
        """
        if(len(readSet)==1): return readSet[0]
        d = dict()
        for r in readSet:
            d[r] = d.get(r,0)+1
        # reverse,sort,reverse
        items = map(list,d.items())
        if len(items)==1: return items[0][0]
        for i in items: i.reverse()
        items.sort(reverse=True)
        for i in items: i.reverse()
        # get the mode and count
        mode = items[0][0]
        if(items[0][1] > items[1][1]): return mode
        return False

    def fitMutProb(self):
        """
        fit a logistic model for the noise/no noise decision
        based on the repeat periods
        """
        # read the relevant data
        data = [(x[0]==x[2],x[0]) for x in self.trainingData if x[0] < TRAIN_MAX_LEN]
        Y = array([int(x[0]) for x in data]) # array of 1/0 if mutated or not
        X = array([x[1] for x in data]) # array of repeat periods
        res = sm.Logit(Y, sm.add_constant(X, prepend=True)).fit()
        self.mutIntercept = res.params[0]
        self.mutSlope = res.params[1]

        # for periods, don't need truncation?
        # fit trunc1 - a contant probability of mutatation for reads >= TRAIN_MAX_LEN
        #data = [int(x[0]==x[2]) for x in self.trainingData if x[0] >= TRAIN_MAX_LEN]
        
        #if(len(data) > MINIMAL_OBS_FOR_TRUNC):
        #    self.trunc1 = sum(data)*1./len(data)
        #else:
        #    data = [int(x[0]==x[2]) for x in self.trainingData]
        #    self.trunc1 = sum(data)*1./len(data)
            
    def fitDecProb(self):
        """
        fit a logistic model for the up/down noise decision
        """
        data = [(x[0]>x[2],x[0]) for x in self.trainingData if x[0] < TRAIN_MAX_LEN\
                and x[0]!=x[2]]
        Y = array([int(x[0]) for x in data])
        X = array([x[1] for x in data])
        res = sm.Logit(Y, sm.add_constant(X, prepend=True)).fit()
        self.decIntercept = res.params[0]
        self.decSlope = res.params[1]
        # fit trunc2 - a contant probability of decrement for reads >= TRAIN_MAX_LEN
        data = [int(x[0]>x[2]) for x in self.trainingData if x[0] >= TRAIN_MAX_LEN\
                and x[0]!=x[2]]
        if(len(data) > MINIMAL_OBS_FOR_TRUNC):
            self.trunc2 = sum(data)*1./len(data)
        else:
            data = [int(x[0]>x[2]) for x in self.trainingData if x[0]!=x[2]]
            self.trunc2 = sum(data)*1./len(data)

    def fitStepProb():
        stepCounts = [10,10,10] # add pseudocounts (1, 2, 3+)
        data = [(abs(x[0]-x[2])-1) for x in self.trainingData if x[0] < TRAIN_MAX_LEN\
            and x[0]!=x[2]]
        for item in data:
            if item > 2: 
                stepCounts[2] += 1
            else:
                stepCounts[item] += 1
        stepCounts = stepCounts/sum(stepCounts)
        self.step1 = stepCounts[0]
        self.step2 = stepCounts[1]
        self.step3 = stepCounts[2]

    def fitPois(self):
        """
        fit a Poisson model for the number of noise steps
        """
        data = [(abs(x[0]-x[2])-1,x[0]) for x in self.trainingData if x[0] < TRAIN_MAX_LEN\
                and x[0]!=x[2]]
        Y = array([int(x[0]) for x in data]) # diff from ref
        X = array([x[1] for x in data]) # STR unit length
        res = sm.Poisson(Y, sm.add_constant(X, prepend=True)).fit()
        self.poisIntercept = res.params[0]
        self.poisSlope = res.params[1]

        # fit trunc3 - a contant mean of mutation for reads >= TRAIN_MAX_LEN
        #data = [abs(x[0]-x[2])-1 for x in self.trainingData if x[0] >= TRAIN_MAX_LEN\
        #        and x[0]!=x[2]]
        #if(len(data) > MINIMAL_OBS_FOR_TRUNC):
        #    self.trunc3 = sum(data)*1./len(data)
        #else:
        #    data = [abs(x[0]-x[2])-1 for x in self.trainingData if x[0]!=x[2]]
        #    self.trunc3 = sum(data)*1./len(data)

    def getTransProb(self,a,b,L):
        """
        given the noise model, what is the probability of observing
        STR=b when the true value is a?
        """
        mutProb = invLogit(self.mutIntercept + self.mutSlope * L)
        #decProb = invLogit(self.decIntercept + self.decSlope * a) 
        poisMean= math.exp(self.poisIntercept + self.poisSlope * L)
        # use the truncation values
        #mutProb = max(mutProb,self.trunc1)
        #decProb = min(decProb,self.trunc2)
        #poisMean = min(poisMean,self.trunc3)
        if(a==b): return mutProb
        if(a!=b):
            diff = abs(a-b)
            #if(a>b): # decrement
                #return (1-mutProb) * decProb * ppois(diff-1,poisMean)
            return (1-mutProb)*ppois(diff-1,poisMean)
            #else: # increment
                #return (1- mutProb) * (1-decProb) * ppois(diff-1,poisMean)
                 
        
    def ReadFromFile(self,filename):
        """
        read the noise model params from a file
        """
        # parse each file and turn into a dictionary
        f = [line.split("=") for line in file(filename).readlines()]
        d = dict([ (line[0],float(line[1])) for line in f])
        self.mutIntercept = d["mutIntercept"]
        self.mutSlope     = d["mutSlope"]
        #self.decIntercept = d["decIntercept"]
        #self.decSlope = d["decSlope"]
        self.poisIntercept = d["poisIntercept"]
        self.poisSlope = d["poisSlope"]
        #self.trunc1 = d["trunc1"]
        #self.trunc2 = d["trunc2"]
        #self.trunc3 = d["trunc3"]
        

    def WriteToFile(self,filename):
        """
        write the noise model parameters to a file
        """
        skel = """mutIntercept=%s
mutSlope=%s
poisIntercept=%s
poisSlope=%s
"""
        f = skel%(self.mutIntercept, self.mutSlope, self.poisIntercept,\
                      self.poisSlope)
        file(filename,"wb").write(f)
#        skel = """mutIntercept=%s
#mutSlope=%s
#decIntercept=%s
#decSlope=%s
#poisIntercept=%s
#poisSlope=%s
#trunc1=%s
#trunc2=%s
#trunc3=%s
#"""
#        f = skel%(self.mutIntercept,self.mutSlope,self.decIntercept,\
#                  self.decSlope,self.poisIntercept,self.poisSlope,\
#                  self.trunc1,self.trunc2,self.trunc3)
#        file(filename,"wb").write(f)

class Genotyper:
    def __init__(self,noiseModel,sex):
        self.noiseModel = noiseModel
        self.currMatSize = MAX_STR_LEN
        self.prepareTransMat()
        self.sex = sex

    def prepareTransMat(self):
        """
        prepare the transition matrix based on the given noise model
        """
        self.transMat = [[[math.log(SMALL_CONST + self.noiseModel.getTransProb(i,j,L))\
                          for j in xrange(self.currMatSize)] for i in xrange(self.currMatSize)] for L in range(2,8)]
        
    def calcLogLik(self,STRa,STRb,reads,L):
        """
        calculate the likelihood of the (STRa,STRb) genotype given the reads
        and the noise model and the STR period
        """
        #probsA = [self.noiseModel.getTransProb(STRa,k,L) for k in range(500)]
        probsA = self.transMat[L-2][STRa]
        #probsB = [self.noiseModel.getTransProb(STRb,k,L) for k in range(500)]
        probsB = self.transMat[L-2][STRb]
        joint = map(lambda x,y: math.log((math.exp(x)+math.exp(y))/2) , probsA,probsB)
        logLik = sum([joint[read] for read in reads])
        return(logLik)

    def findMLE(self,reads,L):
        """
        find the MLE for the given set of reads
        """
        possible = list(set(reads))
        if(len(possible)==1):
            try:
                score = self.calcLogLik(possible[0],possible[0],reads,L)
            except: score = 0
            return([[possible[0],possible[0]],score])
        currBest = [possible[0],possible[0]]
        try:
            currBestScore = self.calcLogLik(possible[0],possible[0],reads,L)
        except: currBestScore = 0
        for i in range(len(possible)):
            for j in range(i+1):
                candidA = possible[i]
                candidB = possible[j]
                try:
                    currScore = self.calcLogLik(candidA,candidB,reads,L)
                except: currScore = 0
                if(currScore > currBestScore):
                    # check if too few minor allele reads
                    # for the heterozygote
                    if candidA != candidB:
                        if reads.count(candidA) > reads.count(candidB):
                            maf = reads.count(candidB)*1.0/len(reads)
                        else:
                            maf = reads.count(candidA)*1.0/len(reads)
                        if maf < MIN_HET_FREQ: continue
                    currBestScore = currScore
                    currBest = [candidA,candidB]
        currBest.sort()
        return [currBest,currBestScore]
    

    def findMLE_single(self,reads,L):
        """
        find the MLE for the given set of reads, assuming there's
        only one STR
        """
        possible = list(set(reads))
        if(len(possible)==1):
            try:
                score = self.calcLogLik(possible[0],possible[0],reads,L)
            except: score = 0
            return([[possible[0]],score])
        currBest = [possible[0],possible[0]]
        try:
            currBestScore = self.calcLogLik(possible[0],possible[0],reads,L)
        except: currBestScore = 0
        for i in range(len(possible)):
            candidA = possible[i]
            try:
                currScore = self.calcLogLik(candidA,candidA,reads,L)
            except:
                currScore = 0
            if(currScore > currBestScore):
                currBestScore = currScore
                currBest = [candidA]
        currBest.sort()
        return [currBest,currBestScore]

    def genotype(self,read_container,output_file):
        """
        genotype the whole file. For each str locus 
        determine the genotype. Then write to a new file.
        """
        f = open(output_file, "w")
        isMale = (self.sex == "M")
        locus = 0
        for str_locus in read_container.aligned_str_map_:
            currReads = GetCopyNumbers(read_container.aligned_str_map_[str_locus])
            currReads.sort()
            resid = currReads[0] - int(currReads[0])
            # floor the reads 
            currReads = map(int,currReads)
            # check that the current transition matrix is large enough
#            if max(currReads) > self.currMatSize:
#                # if not - enlarge it 
#                self.currMatSize = max(currReads)
#                self.prepareTransMat()
          
            currChr = str_locus[0]
            start = str_locus[1]
            stop = str_locus[2]
            repeat = read_container.aligned_str_map_[str_locus][0].opt("XR")
            ref = int(read_container.aligned_str_map_[str_locus][0].opt("XC"))
            if isMale and currChr in [CHR_Y,CHR_X]:
                genotypes  = self.findMLE_single(currReads,len(repeat))
            else:
                try:
                    genotypes = self.findMLE(currReads,len(repeat))
                except: genotypes = [],0
            #genotypes = [x+resid for x in genotypes]
            score = genotypes[1]
            genotypes = genotypes[0]
            if len(genotypes) > 0:
                conflicting = len([item for item in currReads if item not in genotypes])
                not_conflicting = len(currReads) - conflicting
                f.write("\t".join(map(str,[currChr, start, stop, repeat, len(repeat), ref, ",".join(map(str, genotypes)), len(currReads), not_conflicting, conflicting, "/".join(map(str,currReads)), round(score,2) ])) + "\n")


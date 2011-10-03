###################################################
#### This file implements:
#### a) Estimation of the STR noise model
#### b) Required statistical functions
#### c) Genotyping given a noise model
###################################################
## who to hate: David Golan (golandavid@gmail.com)
###################################################
## Note: This is more-or-less a translation from R
## please do not expect proper OOD
###################################################

### required imports
import math,csv
import scikits.statsmodels.api as sm
from numpy import array
 
# used to prevent problems coverting to log scale 
SMALL_CONST = 1e-100
# used as a guess for the maximal length. not critical
MAX_STR_LEN = 100
# other consts
CHR_Y = "chrY"
CHR_X = "chrX"
TRAIN_MAX_LEN = 25 
MINIMAL_OBS_FOR_TRUNC = 10
#### stat utils
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


#### a class for training the model and estimation
class noiseModel(object):
    """
    noiseModel: fit a noise model from a list of files or read the parameters
    from a file. 
    """

    def __init__(self,tabFilenames=None):
        """
        first read the tab files and train the model or
        read the parameters from a file.
        """
        # TODO: fix this to read the file once
        # TODO: add support for multiple files
        # NOTE: the current parsing assumes the current structure of the
        # lobSTR output files.
        if(tabFilenames!=None):
            self.trainFromFiles(tabFilenames)
        else:
            print "initialized without files, please use trainFromFiles or readFromFile"
            
        
            
    def trainFromFiles(self,tabFilenames):
        """
        """
        # we need to train the model
        self.chrs = []
        self.unitSize = []
        self.reads = []
        for filename in tabFilenames:
            print "reading file %s"%filename
            self.readTabFile(filename)
        self.train()        
        
    def readTabFile(self,filename):
        """
        readTabFile reads the STRs from one tab file
        along with all the necessary information
        """
        # TODO: fix this so it reads the file only once
        tbl = csv.reader(file(filename,'rb'),delimiter="\t")
        self.chrs.extend([row[0] for row in tbl])
        tbl = csv.reader(file(filename,'rb'),delimiter="\t")
        self.unitSize.extend([int(row[4]) for row in tbl])
        tbl = csv.reader(file(filename,'rb'),delimiter="\t")
        self.reads.extend([map(float,row[-1].split("/")) for row in tbl])
        
    def hasUniqueMode(self,readSet):
        "gets a list of reads and returns the mode if it is unique and F otherwise"
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

    def train(self):
        """
        estimate the parameters of the noise model
        """
        # get the reads from chrX and chrY with a unique mode
        XY_indices = [x for x in range(len(self.chrs)) if self.chrs[x] in [CHR_Y,CHR_X]\
                      and self.hasUniqueMode(self.reads[x])!= False]
        chrD = {}
        for i in XY_indices: chrD[self.chrs[i]] = chrD.get(self.chrs[i],0)+1
        
        # create data for training
        trainingData = []
        modeD = {}
        for XY_index in XY_indices:
            mode = math.floor(self.hasUniqueMode(self.reads[XY_index]))
            for read in self.reads[XY_index]:
                trainingData.append([mode,self.unitSize[XY_index],math.floor(read)])
        self.trainingData = trainingData
        self.fitMutProb()
        self.fitDecProb()
        self.fitPois()

    def fitMutProb(self):
        """
        fit a logistic model for the noise/no noise decision
        """
        # read the relevant data
        data = [(x[0]==x[2],x[0]) for x in self.trainingData if x[0] < TRAIN_MAX_LEN]

        Y = array([int(x[0]) for x in data])
        X = array([x[1] for x in data])
        res = sm.Logit(Y, sm.add_constant(X, prepend=True)).fit()
        self.mutIntercept = res.params[0]
        self.mutSlope = res.params[1]

        # fit trunc1 - a contant probability of mutatation for reads >= TRAIN_MAX_LEN
        data = [int(x[0]==x[2]) for x in self.trainingData if x[0] >= TRAIN_MAX_LEN]

        if(len(data) > MINIMAL_OBS_FOR_TRUNC):
            self.trunc1 = sum(data)*1./len(data)
        else:
            data = [int(x[0]==x[2]) for x in self.trainingData]
            self.trunc1 = sum(data)*1./len(data)
            
        
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

        
    def fitPois(self):
        """
        fit a Poisson model for the number of noise steps
        """
        data = [(abs(x[0]-x[2])-1,x[0]) for x in self.trainingData if x[0] < TRAIN_MAX_LEN\
                and x[0]!=x[2]]
        Y = array([int(x[0]) for x in data])
        X = array([x[1] for x in data])
        res = sm.Poisson(Y, sm.add_constant(X, prepend=True)).fit()
        self.poisIntercept = res.params[0]
        self.poisSlope = res.params[1]

        # fit trunc3 - a contant mean of mutation for reads >= TRAIN_MAX_LEN
        data = [abs(x[0]-x[2])-1 for x in self.trainingData if x[0] >= TRAIN_MAX_LEN\
                and x[0]!=x[2]]
        print data
        if(len(data) > MINIMAL_OBS_FOR_TRUNC):
            self.trunc3 = sum(data)*1./len(data)
        else:
            data = [abs(x[0]-x[2])-1 for x in self.trainingData if x[0]!=x[2]]
            self.trunc3 = sum(data)*1./len(data)

    def getTransProb(self,a,b):
        """
        given the noise model, what is the probability of observing
        STR=b when the true value is a?
        """
        mutProb = invLogit(self.mutIntercept + self.mutSlope * a)
        decProb = invLogit(self.decIntercept + self.decSlope * a) 
        poisMean= math.exp(self.poisIntercept + self.poisSlope * a)
        # use the truncation values
        mutProb = max(mutProb,self.trunc1)
        decProb = min(decProb,self.trunc2)
        poisMean = min(poisMean,self.trunc3)
        if(a==b): return mutProb
        if(a!=b):
            diff = abs(a-b)
            if(a>b): # decrement
                return (1-mutProb) * decProb * ppois(diff-1,poisMean)
            else: # increment
                return (1- mutProb) * (1-decProb) * ppois(diff-1,poisMean)

        
        
    def readFromFile(self,filename):
        """
        read the noise model params from a file
        """
        # parse each file and turn into a dictionary
        f = [line.split("=") for line in file(filename).readlines()]
        d = dict([ (line[0],float(line[1])) for line in f])
        self.mutIntercept = d["mutIntercept"]
        self.mutSlope     = d["mutSlope"]
        self.decIntercept = d["decIntercept"]
        self.decSlope = d["decSlope"]
        self.poisIntercept = d["poisIntercept"]
        self.poisSlope = d["poisSlope"]
        self.trunc1 = d["trunc1"]
        self.trunc2 = d["trunc2"]
        self.trunc3 = d["trunc3"]
        

    def writeToFile(self,filename):
        """
        write the noise model parameters to a file
        """
        skel = """mutIntercept=%s
mutSlope=%s
decIntercept=%s
decSlope=%s
poisIntercept=%s
poisSlope=%s
trunc1=%s
trunc2=%s
trunc3=%s
"""
        f = skel%(self.mutIntercept,self.mutSlope,self.decIntercept,\
                  self.decSlope,self.poisIntercept,self.poisSlope,\
                  self.trunc1,self.trunc2,self.trunc3)
        file(filename,"wb").write(f)


class genotyper(object):
    def __init__(self,noiseModel):
        """
        """
        self.noiseModel = noiseModel
        self.currMatSize = MAX_STR_LEN
        self.prepareTransMat()

    def prepareTransMat(self):
        """
        prepare the transition matrix based on the given noise model
        """

        self.transMat = [[math.log(SMALL_CONST + self.noiseModel.getTransProb(i,j))\
                          for j in xrange(self.currMatSize)] for i in xrange(self.currMatSize)]
                         
        
    def calcLogLik(self,STRa,STRb,reads):
        """
        calculate the likelihood of the (STRa,STRb) genotype given the reads
        and the noise model
        """
        probsA = self.transMat[STRa]
        probsB = self.transMat[STRb]
        joint = map(lambda x,y: math.log((math.exp(x)+math.exp(y))/2) , probsA,probsB)
        logLik = sum([joint[read] for read in reads])
        return(logLik)

    def findMLE(self,reads):
        """
        find the MLE for the given set of reads
        """
        possible = list(set(reads))
        if(len(possible)==1): return([possible[0],possible[0]])
        currBest = [possible[0],possible[0]]
        currBestScore = self.calcLogLik(possible[0],possible[0],reads)
        for i in range(len(possible)):
            for j in range(i+1):
                candidA = possible[i]
                candidB = possible[j]
                currScore = self.calcLogLik(candidA,candidB,reads)
                if(currScore > currBestScore):
                    currBestScore = currScore
                    currBest = [candidA,candidB]
        currBest.sort()
        return(currBest)
    

    def findMLE_single(self,reads):
        """
        find the MLE for the given set of reads, assuming there's
        only one STR
        """
        possible = list(set(reads))
        if(len(possible)==1): return([possible[0]])
        currBest = [possible[0],possible[0]]
        currBestScore = self.calcLogLik(possible[0],possible[0],reads)
        for i in range(len(possible)):
            candidA = possible[i]
            currScore = self.calcLogLik(candidA,candidA,reads)
            if(currScore > currBestScore):
                currBestScore = currScore
                currBest = [candidA]
        currBest.sort()
        return(currBest)

    def isMale(self,filename):
        """
        read the file and look for chrY entries. Not the most efficient way.
        """
        oldFile = csv.reader(file(filename,'rb'),delimiter="\t")
        chrs = [row[0] for row in oldFile]
        return CHR_X in chrs
    
    def genotype(self,fileForGenotyping,outputFile):
        """
        genotype the whole file - go line by line and add the genotype
        at the end of the line. Then write to a new file.
        """
        newFile = []
        oldFile = csv.reader(file(fileForGenotyping,'rb'),delimiter="\t")
        isMale = self.isMale(fileForGenotyping)
        for row in oldFile:
            currReads = [float(x) for x in row[-1].split("/")]
            resid = currReads[0] - int(currReads[0])
            # floor the reads
            currReads = map(int,currReads)
            # check that the current transition matrix is large enough
            if max(currReads) > self.currMatSize:
                # if not - enlarge it 
                self.currMatSize = max(currReads)
                self.prepareTransMat()
          
            currChr = row[0]
            if isMale and currChr in [CHR_Y,CHR_X]:
                genotypes  = self.findMLE_single(currReads)
            else:
                genotypes = self.findMLE(currReads)
            genotypes = [x+resid for x in genotypes]
            newLine   =  "\t".join(row) + "\t" + ",".join(map(str,genotypes))
            newFile.append(newLine)
        print("saving file...")
        file(outputFile,"wb").write('\n'.join(newFile))
        
                
    
        
    
        
    

    
    
    
    





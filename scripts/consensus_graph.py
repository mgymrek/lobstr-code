# TODO
# if adjusted sequence the same, merge

debug = True
clustalw_exe = "clustalw2"
MATCHSCORE = 1
MISMATCHSCORE = -4
GAPSCORE = -2
EXTENDSCORE = -1
EXTENDREF = 15
MAXDIST = 0.97
K = 2
colors = ["green","blue"]
graph = True
import sys
import os
import random
from Bio.Align.Applications import ClustalwCommandline
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eig
from scipy.cluster.vq import kmeans2
from sklearn.cluster import SpectralClustering

class Consensus:
    """
    Consensus: given the reads aligned to an STR locus, determine
    the consensus sequence of the two alleles present
    """

    def __init__(self, all_reads, ref,name):
        # keep track of the reference allele sequence
        self.ref_allele = ref
        # list of reads
        self.reads = all_reads
        self.name = name
        count = 0
        # get pairwise distances
        self.sim_graph, self.W = self.GetSimilarityGraph()

        # remove errors
        # find errors
        conncomps = nx.connected_component_subgraphs(self.sim_graph)
        nseqs = len(self.reads)
        nodes_to_keep = []
        node_colors = []
        keep = 0
        self.idx = []
        for i in range(len(conncomps)):
            if len(conncomps[i]) >= 0.1*nseqs:
                self.idx.extend([(k,keep) for k in conncomps[i]])
                keep+= 1
        self.idx.sort()
        nodes_to_keep = [item[0] for item in self.idx]
        self.idx = [item[1] for item in self.idx]
        # remove from W
        print "keeping %s components"%keep
        self.W =  self.W[nodes_to_keep,:][:,nodes_to_keep].copy()

        # remake graph and sim matrix
        # TODO try k = 1, k = 2 and get best
        # use SpectralClustering from scikits
        #sc_model = SpectralClustering(k=2).fit(self.W)
        #self.idx = sc_model.labels_
        node_colors = []
        print nodes_to_keep
        for i in range(len(self.reads)):
            if i in nodes_to_keep:
                update_index = nodes_to_keep.index(i)
                try:
                    node_colors.append(colors[self.idx[update_index]])
                except:
                    node_colors.append("gray")
            else:
                node_colors.append("red")

        # visualize this graph
        if graph:
            pos=nx.spring_layout(self.sim_graph)
            edge_labels=dict([((u,v,),round(d['weight'],2)) 
                              for u,v,d in self.sim_graph.edges(data=True)]) 
            plt.clf()
            nx.draw(self.sim_graph,pos, node_color = node_colors)
            nx.draw_networkx_edge_labels(self.sim_graph,pos,edge_labels=edge_labels) 
            plt.savefig("/var/www/lobSTR/test/graphtest_%s.png"%self.name)

        self.reads = [self.reads[i] for i in nodes_to_keep]
        if debug:
            for i in range(len(self.reads)): 
                print i, self.idx[i], self.reads[i]

        # Align reads in each group
        for k in range(K):
            readsk = [self.reads[j] for j in range(len(self.reads)) if self.idx[j] == k]
            alignmentk = self.AlignReads(readsk)
            if ref != "":
                seqk_adjusted, seqk_cigar = self.AdjustAlignment(alignmentk)
            else:
                seqk_adjusted,seqk_cigar = "NA","NA"
            if debug:
                print
                print "K:",k
                print "alignment", alignmentk
                print "adjusted", seqk_adjusted
                print "ref     ", self.ref_allele
                print "cigar", seqk_cigar
        

    def GetSimilarity(self,seq1,seq2):
        """
        Get alignment distance between two sequences
        Use Clustal with no gap ext or end gap penalty
        """
        # Clustalw alignment of ref vs. allele
        randint = random.randint(0,10000)
        f = open("%s.fa"%randint,"w")
        f.write(">seq1\n")
        f.write(seq1+"\n")
        f.write(">seq2\n")
        f.write(seq2+"\n")
        f.close()
        clustalw_cline = ClustalwCommandline(clustalw_exe,infile="%s.fa"%randint,endgaps=True,gapopen=8,pairgap=8,pwgapopen=8,gapext=0.05)
        clustalw_cline()
        seq1align, seq2align = self.ParseClustalGetSeqs("%s.aln"%randint)
        os.system("rm %s*"%randint)
        score = self.GetAlignmentScore(seq1align, seq2align)
        return score

    def GetAlignmentScore(self,s1, s2):
        """
        Get score for pairwise alignment
        """
        score = 0
        i = 0
        # get rid of end gaps
        while s1[i] == "-" or s2[i] == "-":
            i += 1
        startind = i
        i = len(s1)-1
        while s1[i] == "-" or s2[i] == "-":
            i -= 1
        endind = i
        overlap = endind-startind+1
        last = "none"
        for i in range(startind,endind+1):
            if s1[i] == s2[i]:
                score += MATCHSCORE
                last = "match"
            elif s1[i] == "-":
                if last == "g1":
                    score += EXTENDSCORE
                else: score += GAPSCORE
                last = "g1"
            elif s2[i] == "-":
                if last == "g2":
                    score += EXTENDSCORE
                else: score += GAPSCORE
                last = "g2"
            elif s1[i] == "N" or s2[i] == "N":
                last = "match"
                overlap -= 1
            else:
                score += MISMATCHSCORE
                last = "mismatch"
        return score*1.0/overlap
            
    def GetSimilarityGraph(self):
        """
        Given a set of reads, return pairwise distances 
        based on alignments
        """
        distgraph = nx.Graph()
        nseqs = len(self.reads)
        W = np.zeros((nseqs,nseqs))
        for i in range(nseqs):
            distgraph.add_node(i)
            for j in range(i,nseqs):
                sim = self.GetSimilarity(self.reads[i],self.reads[j])
                if sim > MAXDIST:
                    if j != i:
                        distgraph.add_edge(i,j,weight=sim)
                        W[i][j] = 1
                        W[j][i] = 1
        return distgraph, W
    
    def ParseClustalGetSeqs(self, alnfile):
        """ 
        Parse clustal .aln file to get
        ref seq and allele seq
        """
        refseq = ""
        alleleseq = ""
        f = open(alnfile,"r")
        line = f.readline()
        block_seqs = []
        while line != "":
            if "CLUSTAL" in line:
                line = f.readline()
                while line.strip() == "":
                    line = f.readline()
            if line != "" and line.strip() == "":
                # process last block and start new one
                if len(block_seqs)>1:
                    refseq += block_seqs[0]
                    alleleseq += block_seqs[1]
                block_seqs = []
            else:
                if "*" not in line:
                    block_seqs.append(line.strip().split()[1])
            line = f.readline()
        # parse last block
        if len(block_seqs) > 1:
            refseq += block_seqs[0]
            alleleseq += block_seqs[1]
        return alleleseq, refseq

    def AdjustAlignment(self, sequence):
        """
        Adjust alignment to have the same start and end coordinates as 
        the reference allele.
        """
        # Clustalw alignment of ref vs. allele
        randint = random.randint(0,10000)
        f = open("%s.fa"%randint,"w")
        f.write(">ref\n")
        f.write(self.ref_allele+"\n")
        f.write(">allele\n")
        f.write(sequence+"\n")
        f.close()
        clustalw_cline = ClustalwCommandline(clustalw_exe,infile="%s.fa"%randint,endgaps=True,gapopen=8,pairgap=5,pwgapopen=8,gapext=0.0)
        clustalw_cline()
        allele, ref = self.ParseClustalGetSeqs("%s.aln"%randint)
        os.system("rm %s*"%randint)
        if debug:
            print "allele",allele
            print "ref   ",ref
        # remove end gaps from ref
        i = 0
        while ref[i] == "-":
            i = i + 1
        ref = ref[i:]
        allele = allele[i:]
        i = len(allele)-1
        while ref[i] == "-":
            i = i-1
        ref = ref[:i+1]
        allele = allele[:i+1]
        # remove extra 10 flanking added to ref
        numrev = 0
        i = 0
        while numrev < EXTENDREF:
            if ref[i] != "-":
                numrev = numrev + 1
            i = i+1
        ref = ref[i:]
        allele = allele[i:]
        i = len(ref)-1
        numrev = 0
        while numrev < EXTENDREF:
            if ref[i] != "-":
                numrev = numrev + 1
            i = i-1
        ref = ref[:i+1]
        allele = allele[:i+1]

        # make cigar score (use m to indicate SNP)
        cigar = []
        for i in range(len(allele)):
            if allele[i] == ref[i]: cigar.append("M")
            elif allele[i] == "-": cigar.append("D")
            elif ref[i] == "-": cigar.append("I")
            else: cigar.append("m")
        if len(cigar) == 0: return "",""
        if len(cigar) == 1: return allele, cigar[0]
        cigarstring = ""
        curcount = 1
        curchar = cigar[0]
        for i in range(1,len(cigar)):
            newchar = cigar[i]
            if newchar != curchar:
                cigarstring = cigarstring + str(curcount)+curchar
                curcount = 1
                curchar = newchar
            else:
                curcount = curcount + 1
        cigarstring = cigarstring + str(curcount)+curchar
        allele = allele.replace("-","")
        return allele, cigarstring

    def AlignReads(self, read_list):
        """
        Perform multiple alignment on a list of sequences
        """
        if len(read_list) == 1: return read_list[0]
        if len(read_list) == 0: return ""
        # make fasta file of sequences
        randint = random.randint(0,1000)
        fname="%s.fa"%randint
        f = open(fname,"w")
        for i in range(len(read_list)):
            f.write(">seq%s\n"%i)
            f.write(read_list[i]+"\n")
        f.close()
        clustalw_cline = ClustalwCommandline(clustalw_exe,infile="%s.fa"%randint)
        clustalw_cline()
        alignment = self.ParseClustal("%s.aln"%randint)
        os.system("rm %s*"%randint)
        return alignment

    def ParseClustal(self, alnfile):
        """ 
        Parse clustal .aln file to get
        consensus sequence
        """
        seq = ""
#        if debug: print(open(alnfile,"r").read())
        f = open(alnfile,"r")
        line = f.readline()
        block_seqs = []
        while line != "":
            if "CLUSTAL" in line:
                line = f.readline()
                while line.strip() == "":
                    line = f.readline()
            if line != "" and line.strip() == "":
                # process last block and start new one
                if len(block_seqs)>0:
                    for i in range(len(block_seqs[0])):
                        chars = []
                        for b in range(len(block_seqs)):
                            newchar = block_seqs[b][i]
                            chars.append(newchar)
                        conchar = self.GetMajorityVote(chars)
                        if conchar != "-":
                            seq += conchar
                block_seqs = []
            else:
                if "*" not in line:
                    block_seqs.append(line.strip().split()[1])
            line = f.readline()
        # parse last block
        if len(block_seqs) > 0:
            for i in range(len(block_seqs[0])):
                chars = []
                for b in range(len(block_seqs)):
                    newchar = block_seqs[b][i]
                    chars.append(newchar)
                conchar = self.GetMajorityVote(chars)
                if conchar != "-":
                    seq += conchar
        return seq
    
    def GetMajorityVote(self, charlist):
        """
        Get majority vote from a list of characters
        """
        if len(charlist) == 0: return ''
        maxitem = ''
        maxcount = 0
        for item in set(charlist):
            num = charlist.count(item)
            if num>=maxcount and item !="-":
                maxitem = item
                maxcount = num
        return maxitem
        
 

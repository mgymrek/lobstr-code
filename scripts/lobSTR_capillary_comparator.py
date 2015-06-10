#!/usr/bin/env python
"""
Compare capillary vs. lobSTR calls

This script is part of lobSTR_validation_suite.sh and
is not mean to be called directly.
"""

import argparse
import numpy as np
import pandas as pd
import sys
from scipy.stats import pearsonr

def ConvertSample(x):
    """
    Convert HGDP samples numbers to standard format
    HGDPXXXXX
    """
    num = x.split("_")[1]
    zeros = 5-len(num)
    return "HGDP"+"0"*zeros + num

def LoadCapillaryFromStru(capfile, convfile):
    """
    Input:
      capfile: filename for .stru file
      convfile: filename giving illumina sample id->HGDP id
    Output: data frame with capillary calls

    Construct data frame with:
       marker
       sample
       allele1.cap
       allele2.cap

    Use sample names converted to Illumina format
    Ignore alleles that are -9,-9
    """
    # Load conversions
    conv = pd.read_csv(convfile, sep="\t")
    converter = dict(zip(conv.hgdp, conv.sample))
    # Load genotypes
    markers = []
    samples = []
    allele1s = []
    allele2s = []
    f = open(capfile, "r")
    marker_names = f.readline().strip().split()
    line = f.readline()
    while line != "":
        items = line.strip().split()
        ident = "HGDP_%s"%items[0]
        pop_code = items[1]
        pop_name = items[2]
        geo = items[3]
        geo2 = items[4]
        alleles1 = items[5:]
        line = f.readline() # get second allele for the individual
        items = line.strip().split()
        ident2 = "HGDP_%s"%items[0] 
        if ident != ident2:
            sys.stderr.write("ERROR parsing .stru file for individual %s\n"%items[0])
            sys.exit(1)
        alleles2 = items[5:]
        sample = converter.get(ConvertSample(ident), "NA")
        if sample != "NA":
            for i in range(len(alleles1)):
                if str(alleles1[i]) != "-9" and str(alleles2[i]) != "-9":
                    markers.append(marker_names[i])
                    samples.append(sample)
                    allele1s.append(int(alleles1[i]))
                    allele2s.append(int(alleles2[i]))
        line = f.readline()
    return pd.DataFrame({"marker": markers, "sample": samples, \
                         "allele1.cap": allele1s, "allele2.cap": allele2s})

def GetAllele(x, allele_num):
    if allele_num == 1:
        al = x["allele1.cap"]
    else: al = x["allele2.cap"]
    raw_allele = ((al-x["effective_product_size"])/x["period"])*x["period"]
    corr_allele = raw_allele - x["correction"]
    return corr_allele

def GetDosage(a1, a2):
    if str(a1) == "." or str(a2) == ".": return "NA"
    else: return (float(a1)+float(a2))*0.5

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lobSTR", help="Tab file with lobSTR calls created by lobSTR_vcf_to_tab.py", type=str, required=True)
    parser.add_argument("--cap", help=".stru file with capillary calls", type=str, required=True)
    parser.add_argument("--corrections", help="Tab file with Marshfield marker corrections", type=str, required=True)
    parser.add_argument("--sample-conversions", help="Tab file with conversion between sample ids", type=str, required=True)
    parser.add_argument("--output-stats", help="Output sample, call, and locus level stats to files with this prefix", type=str, required=False)
    args = parser.parse_args()

    LOBFILE = args.lobSTR
    CAPFILE = args.cap
    CORRFILE = args.corrections
    CONVFILE = args.sample_conversions

    # Load lobSTR calls and corrections
    lob = pd.read_csv(LOBFILE, sep="\t")
    corr = pd.read_csv(CORRFILE, sep="\t")

    # Load capillary calls to data frame
    cap = LoadCapillaryFromStru(CAPFILE, CONVFILE)

    # Merge ddatasets
    res = pd.merge(lob, corr, on=["chrom", "start"])
    res = pd.merge(res, cap, on=["marker", "sample"])
    res["allele1.cap.corr"] = res.apply(lambda x: GetAllele(x, 1), 1)
    res["allele2.cap.corr"] = res.apply(lambda x: GetAllele(x, 2), 1)
    res["correct"] = res.apply(lambda x: str(x["allele1"])==str(x["allele1.cap.corr"]) and \
                               str(x["allele2"])==str(x["allele2.cap.corr"]), 1)
    res["dosage_lob"] = res.apply(lambda x: GetDosage(x["allele1"], x["allele2"]), 1)
    res["dosage_cap"] = res.apply(lambda x: GetDosage(x["allele1.cap.corr"], x["allele2.cap.corr"]), 1)

    ##### Stats #####
    if args.output_stats:
        # Call level stats
        res.to_csv(args.output_stats+".calllevel.tab", index=False, sep="\t")
        # Sample level stats
        sample_level = res.groupby("sample", as_index=False).agg({"DP": np.mean,
                                                                  "Q": np.mean,
                                                                  "start": len,
                                                                  "correct": np.mean})
        sample_level.to_csv(args.output_stats+".samplelevel.tab", index=False, sep="\t")
        # Locus level stats
        res["length"] = res["end_x"]-res["start"]+1
        locus_level = res.groupby(["chrom","start"], as_index=False).agg({"length": np.mean,
                                                                          "DP": np.mean,
                                                                          "Q": np.mean,
                                                                          "GT": len,
                                                                          "SB": np.mean,
                                                                          "DISTENDS": np.mean,
                                                                          "correct": np.mean})
        locus_level.to_csv(args.output_stats+".locuslevel.tab", index=False, sep="\t")

    ##### Results #####
    sys.stdout.write("########## Results ########\n")
    # Get stats about calls
    num_samples = len(set(res[res["allele1"].apply(str)!="."]["sample"]))
    num_markers = len(set(res[res["allele1"].apply(str)!="."]["marker"]))
    num_nocalls = res[res["allele1"].apply(str)=="."].shape[0]
    sys.stdout.write("# Samples: %s\n"%num_samples)
    sys.stdout.write("# Markers: %s\n"%num_markers)
    sys.stdout.write("# No call rate: %s\n"%(num_nocalls*1.0/res.shape[0]))
    sys.stdout.write("# Number of calls compared: %s\n"%res.shape[0])

    # Accuracy
    acc = np.mean(res[res["allele1"].apply(str)!="."]["correct"])
    sys.stdout.write("# Accuracy: %s\n"%acc)

    # R2
    dl = map(float, list(res[res["allele1"].apply(str)!="."]["dosage_lob"]))
    dc = map(float, list(res[res["allele1"].apply(str)!="."]["dosage_cap"]))
    r2 = pearsonr(dl, dc)[0]**2
    sys.stdout.write("# R2: %s\n"%r2)

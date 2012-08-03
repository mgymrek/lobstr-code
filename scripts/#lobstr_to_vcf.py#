def usage():
    print """
python lobstr_to_vcf.py [OPTIONS] --gen <lobstr genotypes.tab file> --sample <STRING> --out <STRING>  [--zip]

REQUIRED:
--gen: lobSTR allelotype .genotypes.tab output
--sample: string denoting the name of this sample
--out: prefix to name output file

OPTIONS:
--zip: file is in gzip format
--exclude: file of "chrom\\tpos" of positions to exclude. 
           For downstream analysis, it is beneficial to exlucde
           any STRs with the same starting point but different motifs,
           as this will cause errors when using vcftools.
--help,-h: print this message
--verbose,-v: print helpful status messages

"""

import os
import sys
import getopt
import vcf
import datetime
import math
import gzip

try:
    opts, args = getopt.getopt(sys.argv[1:], "hv", ["help","gen=", "sample=", "out=","exclude=","zip"])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)
args = [item[0] for item in opts]

if (not("--gen" in args and "--out" in args and "--sample" in args)):
    usage()
    sys.exit(2)

# initialize variables
genfile = ""
output_prefix = ""
sample = ""
verbose = False
gzipped = False
excludefile = ""


for o,a in opts:
    if o == "-v": verbose = True
    if o == "--out": output_prefix = a
    if o == "--gen": genfile = a
    if o == "--sample": sample = a
    if o == "--exclude": excludefile = a
    if o == "--zip": gzipped = True
    if o == "-h" or o == "--help":
        usage()
        sys.exit(0)

###############################
# methods
def GetVCFHeader(genfile, output_prefix):
    """ Get the VCF header string """
    if gzipped:
        header_lines = gzip.open(genfile,"rb").readlines()[0:2]
    else:
        header_lines = open(genfile,"r").readlines()[0:2]        
    ref = ""
    if "index-prefix=" in header_lines[0]:
        hline = header_lines[0].split(";")[:-1]
        newhline = []
        for item in hline:
            if "=" not in item:
                newhline.append(item+"=")
            else: newhline.append(item)
        alignment_objects = dict([(item.split("=")[0],item.split("=")[1]) for item in newhline])
        ref = alignment_objects["index-prefix"]
    header = open("%s.header.vcf"%output_prefix,"w")
    header.write("##fileformat=VCFv4.1\n")
    header.write("##fileDate=%s\n"%str(datetime.datetime.now()).split()[0].replace("-",""))
    if "#" in header_lines[0]:
        header.write("##source=%s:%s\n"%(header_lines[0].replace("#","").replace("=",":").strip(),header_lines[1].replace("#","").replace("=",":").strip()))
    if ref != "":
        header.write("##reference=file://%s\n"%ref.replace(".","").strip())
    header.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n")
    header.write("##INFO=<ID=PERIOD,Number=1,Type=Integer,Description=\"Period of the repeat motif\">\n")
    header.write("##INFO=<ID=MOTIF,Number=1,Type=String,Description=\"Repeat motif\">\n")
    header.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
    header.write("##FORMAT=<ID=SUPP,Number=1,Type=Integer,Description=\"Number of reads supporting genotype call\">\n")
    header.write("##FORMAT=<ID=CONFLICT,Number=1,Type=Integer,Description=\"Number of reads conflicting genotype call\">\n")
    header.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    header.write("##FORMAT=<ID=GB,Number=1,Type=String,Description=\"Genotype given in bp difference from reference\">\n")
    header.write("##FORMAT=<ID=ALL,Number=1,Type=String,Description=\"All alleles seen\">\n")
    header.write("##FORMAT=<ID=S1,Number=1,Type=Float,Description=\"Allele 1 score\">\n")
    header.write("##FORMAT=<ID=S2,Number=1,Type=Float,Description=\"Allele 2 score\">\n")
    header.write("##ALT=<ID=STRVAR,Description=\"Short tandem variation\">\n")
    header.write("##INFO=<ID=REF,Number=1,Type=Float,Description=\"Reference copy number\">\n")
    header.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s"%sample.strip())
    header.close()
    return header.name

def ProcessGenotypes(genfile, vcf_writer, pos_to_exclude):
    """ print vcf record for each allelotyped locus"""
    if gzipped:
        f = gzip.open(genfile,"rb")
    else:
        f = open(genfile, "r")
    # get rid of headers
    line = f.readline()
    while "#" in line: line = f.readline()
    while line.strip() != "" and len(line.split("\t")) >= 15:
        genotype_failed = False
        items = line.strip().split("\t")
        chrom = items[0]
        start = int(items[1])
        if (chrom, start) in pos_to_exclude: 
            genotype_failed = True
        end = int(items[2])
        motif = items[3]
        period = int(items[4])
        ref_copy_num = float(items[5])
        if ref_copy_num <= 0: genotype_failed = True
        allelotype = items[6]
        alleles = [allelotype.split(",")[0],allelotype.split(",")[1]]
        coverage = int(items[7])
        agree = int(items[8])
        conflict = int(items[9])
        all_alleles = items[10]
        score = float(items[11])
        score1 = float(items[12])
        score2 = float(items[13])
        try:
            partial_coverage = int(items[12])
            max_partial_string = items[13]
            partial_read_string = items[14]
        except:
            partial_coverage = 0
        CHROM = chrom
        POS = start
        ID = "."
        QUAL = score
        FILTER = ""
        INFO = {"REF":ref_copy_num, "END":end, "PERIOD":len(motif.strip()), "MOTIF": motif.strip()}
        FORMAT="GT:GB:DP:SUPP:CONFLICT:ALL:S1:S2"
        SAMPLE_INDICES = {sample:0}
        genotype_bp = alleles[0]+"/"+alleles[1]
        genotype = ""
        REF = motif[0]
        ALT = dict([("<STRVAR>",0) for item in alleles if (item != "0" and item != "NA")]).keys()
        if len(ALT) == 0: ALT = ["."]
        if alleles[0] == "0":
            genotype += "0/"
        elif alleles[0] == "NA":
            genotype += "./"
        else:
            genotype += "1/"
        if alleles[1] == "0":
            genotype += "0"
        elif alleles[1] == "NA":
            genotype += "."
        else:
            genotype += "1"
        if partial_coverage != 0 and coverage == 0:
            # deal with partial coverage TODO
            genotype_failed = True
        if not genotype_failed:
            all_alleles = all_alleles.replace(":","|")
            call = vcf.parser._Call("",sample, data={"GT":genotype,"GB":genotype_bp,"DP":coverage,"SUPP":agree,"CONFLICT":conflict,"ALL":all_alleles,"S1":score1,"S2":score2})
            SAMPLES = [call]
            if not genotype_failed:
                vcf_record = vcf.parser._Record(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE_INDICES, samples=SAMPLES)
                vcf_writer.write_record(vcf_record)
        line = f.readline()

def GetPositionsToExclude(efile):
    pos = []
    if efile == "": return pos
    f = open(efile,"r")
    line = f.readline()
    while line != "":
        chrom, start = line.strip().split("\t")[0:2]
        pos.append((chrom,int(start)))
        line = f.readline()
    f.close()
    return pos
###############################
# main

def main():
    # prepare positions to exclude
    if verbose: print("Getting positions to exclude...")
    pos_to_exclude = GetPositionsToExclude(excludefile)

    # prepare VCF file with header
    if verbose: print("Setting VCF header...")
    header = GetVCFHeader(genfile, output_prefix)
    vcf_reader = vcf.Reader(open(header,"r"))
    vcf_output_file = open("%s.vcf"%output_prefix,"w")
    vcf_writer = vcf.Writer(vcf_output_file, vcf_reader)

    # populate with genotype information
    if verbose: print("Setting allelotype information...")
    ProcessGenotypes(genfile, vcf_writer, pos_to_exclude)
    vcf_output_file.close()

    # fix trailing white spaces and move fileformat to top
    lines = open("%s.vcf"%output_prefix,"r").readlines()
    f = open("%s.vcf"%output_prefix,"w")
    f.write("##fileformat=VCFv4.1\n")
    for line in lines: 
        if "fileformat" not in line:
            f.write(line.strip()+"\n")
    f.close()
    

main()



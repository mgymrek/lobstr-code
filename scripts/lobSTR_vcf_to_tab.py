#!/usr/bin/env python

import os
import sys
import vcf

try:
    VCFFILE = sys.argv[1]
except:
    sys.stderr.write("Usage: vcf_to_tab.py VCFFILE\n")
    sys.exit(1)

if not os.path.exists(VCFFILE):
    sys.stderr.write("ERROR: %s does not exist\n"%VCFFILE)
    sys.exit(1)

reader = vcf.Reader(open(VCFFILE, "rb"))

sys.stdout.write("\t".join(["chrom","start","sample","allele1","allele2","cov","qual","reads"])+"\n")

for record in reader:
    chrom = record.CHROM
    start = record.POS
    for sample in record:
        name = sample.sample
        if sample["GT"]:
            a1, a2 = sample["GB"].split("/")
            cov = sample["DP"]
            qual = sample["Q"]
            reads = sample["ALLREADS"]
        else:
            a1, a2 = ".", "."
            cov = 0
            qual = -1
            reads = "."
        sys.stdout.write("\t".join(map(str, [chrom, start, name, a1, a2, cov, qual, reads]))+"\n")

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

format_keys = reader.formats.keys()
# Remove keys we don't want displayed as extra tabs
if "GB" in format_keys: format_keys.remove("GB")
if "PL" in format_keys: format_keys.remove("PL")
if "GL" in format_keys: format_keys.remove("GL")
sys.stdout.write("\t".join(["chrom","start","end","sample","allele1","allele2"] + format_keys)+"\n")

for record in reader:
    chrom = record.CHROM
    start = record.POS
    end = record.INFO["END"]
    for sample in record:
        name = sample.sample
        vals = []
        if sample["GT"]:
            vals.extend(sample["GB"].split("/"))
            for fk in format_keys:
                vals.append(sample[fk])
        else:
            vals = [".","."]
            for fk in format_keys:
                vals.append(0)
        sys.stdout.write("\t".join(map(str, [chrom, start, end, name] + vals))+"\n")

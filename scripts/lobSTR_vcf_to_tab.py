#!/usr/bin/env python

import os
import sys
import vcf

try:
    VCFFILE = sys.argv[1]
except:
    sys.stderr.write("Usage: vcf_to_tab.py VCFFILE\n")
    sys.exit(1)

if not os.path.exists(VCFFILE) and VCFFILE != "-":
    sys.stderr.write("ERROR: %s does not exist\n"%VCFFILE)
    sys.exit(1)

if VCFFILE == "-":
    reader = vcf.Reader(sys.stdin)
else:
    reader = vcf.Reader(open(VCFFILE, "rb"))

format_keys = reader.formats.keys()
# Remove keys we don't want displayed as extra tabs
if "GB" in format_keys: format_keys.remove("GB")
if "PL" in format_keys: format_keys.remove("PL")
if "GL" in format_keys: format_keys.remove("GL")
sys.stdout.write("\t".join(["chrom","start","end","sample","allele1","allele2"] + format_keys + ["motif","filter_locus","filter_call"])+"\n")

for record in reader:
    chrom = record.CHROM
    start = record.POS
    end = record.INFO["END"]
    motif = record.INFO["MOTIF"]
    filter_locus = record.FILTER
    if filter_locus is None or len(filter_locus) == 0: filter_locus = ["PASS"]
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
                try:
                    v = sample[fk]
                except:
                    v = 0
                vals.append(v)
        try:
            filter_call = sample["FT"]
        except: filter_call = "PASS"
        sys.stdout.write("\t".join(map(str, [chrom, start, end, name] + vals + [motif, ",".join(filter_locus), filter_call]))+"\n")

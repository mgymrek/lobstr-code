#!/usr/bin/env python

import os
import sys
import vcf

try:
    VCFFILE = sys.argv[1]
    LOCUS = sys.argv[2]
    SAMPLE = sys.argv[3]
except:
    sys.stderr.write("Usage: vcf_extract.py <VCFFILE> <LOCUS> <SAMPLE>\nLocus has format chr:pos.\n")
    sys.exit(1)

if not os.path.exists(VCFFILE):
    sys.stderr.write("ERROR: %s does not exist\n"%VCFFILE)
    sys.exit(1)
if not os.path.exists(VCFFILE + ".tbi"):
    sys.stderr.write("ERROR: %s does not have an index. Please run tabix -p vcf %s\n"%(VCFFILE, VCFFILE))
    sys.exit(1)

try:
    vcfreader = vcf.Reader(open(VCFFILE, "rb"))
except:
    sys.stderr.write("ERROR: Could not open %s. Is this a valid VCF file?\n"%VCFFILE)
    sys.exit(1)

try:
    chrom, start = LOCUS.split(":")
except:
    sys.stderr.write("ERROR: Incorrect locus format. Must take the form chr:pos. e.g. chr1:10000.\n")
    sys.exit(1)
try:
    start = int(start)
except ValueError:
    sys.stderr.write("ERROR: Start coordinate must be an integer\n")
    sys.exit(1)

if SAMPLE not in vcfreader.samples:
    sys.stderr.write("ERROR: Sample %s not found in VCF\n"%SAMPLE)
    sys.exit(1)

try:
    record = vcfreader.fetch(chrom, start)
except IOError:
    sys.stderr.write("ERROR: Error fetching record. Make sure VCF file is indexed (tabix -p vcf <VCFFILE>)\n")
    sys.exit(1)
except ValueError:
    sys.stderr.write("ERROR: Error fetching record. This record likely does not exist\n")
    sys.exit(1)

found_sample = False
for sample in record.samples:
    if sample.sample == SAMPLE:
        found_sample = True
        if sample["GT"] is not None:
            sys.stdout.write("%s %s\n"%(LOCUS, SAMPLE))
            sys.stdout.write("Genotype: %s\n"%sample["GT"])
            try:
                sys.stdout.write("GB: %s\n"%sample["GB"])
            except: pass
            try:
                sys.stdout.write("DP: %s\n"%sample["DP"])
            except: pass
            try:
                sys.stdout.write("Q: %s\n"%sample["Q"])
            except: pass
            try:
                sys.stdout.write("ALLREADS: %s\n"%sample["ALLREADS"])
            except: pass
        else:
            sys.stderr.write("No call\n")
            sys.exit(1)

if not found_sample:
    sys.stderr.write("ERROR: Did not find sample in VCF. This shouldn't happen.\n")
    sys.exit(1)

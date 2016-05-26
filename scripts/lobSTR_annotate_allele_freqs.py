#!/usr/bin/env python
"""
Add a field to VCF with population allele frequencies based 
on a reference VCF
"""

import argparse
import gzip
import math
import numpy as np
import os
import sys
import vcf
import tempfile

def MSG(string):
    sys.stderr.write(string.strip() + "\n")

def GetWriter(reader):
    """
    Get VCF Writer with the appropriate metadata
    """
    tmpdir = tempfile.mkdtemp(prefix="lobstr.")
    tmpfile = os.path.join(tmpdir, "header.vcf")
    f = open(tmpfile, "w")
    for line in reader._header_lines: f.write(line.strip() + "\n")
    f.write("##FORMAT=<ID=AF,Number=1,Type=String,Description=\"Population allele frequencies.\">\n")
    f.write("#" + "\t".join(reader._column_headers + reader.samples) + "\n")
    f.close()
    writer = vcf.Writer(sys.stdout, vcf.Reader(open(tmpfile, "rb")))
    return writer

def GetAfreqs(chrom, pos, refreader):
    # Return dictionary of allele->freq
    try:
        records = refreader.fetch(chrom, pos, pos+1)
    except: return {}
    for record in records:
        if record.CHROM == chrom and record.POS == pos:
            if record.ALT[0] is None: return {0: 1.0}
            allele_freqs = [1-sum(record.aaf)] + record.aaf
            allele_sizes = [0] + [len(record.ALT[i])-len(record.REF) for i in range(len(record.ALT))]
            return dict(zip(allele_sizes, allele_freqs))
    return {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="Input VCF to be annotated", type=str, required=True)
    parser.add_argument("--refpop", help="VCF of reference population. Must be indexed and bgzipped", type=str, required=True)
    args = parser.parse_args()

    # Open ref population
    refreader = vcf.Reader(open(args.refpop, "rb"))

    # Open VCF reader
    if args.vcf == "-":
        try:
            reader = vcf.Reader(sys.stdin)
        except (IOError, ValueError, StopIteration) as e:
            MSG("Problem reading VCF file. Is this a valid VCF?")
            sys.exit(1)
    else:
        try:
            reader = vcf.Reader(open(args.vcf, "rb"))
        except (IOError, ValueError) as e:
            MSG("Problem reading VCF file. Is this a valid VCF?")
            sys.exit(1)

    # Set up VCF writer with new filter FORMAT field
    writer = GetWriter(reader)

    for record in reader:
        # Set up afreq info
        if "AF" not in reader.formats:
            record.add_format("AF")
        samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(":"))
        for fmt in samp_fmt._fields:
            if fmt == "AF":
                entry_type = "String"
                entry_num = 1
            else:
                entry_type = reader.formats[fmt].type
                entry_num = reader.formats[fmt].num
            samp_fmt._types.append(entry_type)
            samp_fmt._nums.append(entry_num)
        new_samples = []
        # Look up afreqs (dictionary of genotype->frequency
        afreqs = GetAfreqs(record.CHROM, record.POS, refreader)
        # Get sample info
        for sample in record:
            if sample["GT"]:
                gt1, gt2 = map(float, sample["GT"].split("/"))
                afreqstring = "%.3f/%.3f"%(afreqs.get(gt1, -1), afreqs.get(gt2, -1))
            else: 
                afreqstring = "-1/-1"
            sampdat = []
            for i in range(len(samp_fmt._fields)):
                key = samp_fmt._fields[i]
                if key != "AF":
                    sampdat.append(sample[key])
                else: sampdat.append(afreqstring)
            call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
            new_samples.append(call)
        record.samples = new_samples
        writer.write_record(record)
        

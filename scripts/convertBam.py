def usage():
    print """
python convertBam.py <str.table> <input.sorted.bam> <output.bam>

converts a standard bam format to include the flags needed for lobSTR genotyping
"""

import os
import sys
import pysam

try:
    str_table = sys.argv[1]
    input_bam = sys.argv[2]
    output_bam = sys.argv[3]
except:
    usage()
    sys.exit(2)


# open bam files
input_bam = pysam.Samfile(input_bam, 'rb')
output_bam = pysam.Samfile(output_bam, 'wb', template = input_bam)

# for each STR
strs = open(str_table, "r")
line = strs.readline()
while (line != ""):
    # - fetch all reads overlapping region
    items = line.split("\t")
    chrom, start, end, repeat, ref = items[0],items[1], items[2], items[15], items[5]
    for aligned_read in input_bam.fetch(chrom, int(start), int(end)):
        # flags
        diff = 0
        for cigar_op in aligned_read.cigar:
            op, length = cigar_op
            if op == 1:
                diff += length
            elif op == 2:
                diff -= length
        if abs(diff) > 10: print aligned_read.cigar
        tags = []
        tags.append(("XS", int(start)))
        tags.append(("XE", int(end)))
        tags.append(("XR", repeat.strip()))
        tags.append(("XC", float(ref)))
        tags.append(("XD", diff))
        aligned_read.tags = tags
        # make sure has at least some flanking region TODO
        output_bam.write(aligned_read)
    line = strs.readline()

strs.close()
input_bam.close()
output_bam.close()

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
    input_bam_file = sys.argv[2]
    output_bam_file = sys.argv[3]
except:
    usage()
    sys.exit(2)


# open bam files
input_bam = pysam.Samfile(input_bam_file, 'rb')
output_bam = pysam.Samfile("temp", 'wb', template = input_bam)

# for each STR
strs = open(str_table, "r")
line = strs.readline()
while (line != ""):
    # - fetch all reads overlapping region
    items = line.split("\t")
    chrom, start, end, repeat, ref = items[0],items[1], items[2], items[14], items[4]
    for aligned_read in input_bam.fetch(chrom, int(start), int(end)):
        # make sure it spans the entire repeat
        if (aligned_read.qstart <= int(end) and aligned_read.qstart >= int(start)) \
                or (aligned_read.qend <=int(end) and aligned_read.qend >=int(start)):
            continue # doesn't span whole STR region
        # flags
        diff = 0
        for cigar_op in aligned_read.cigar:
            op, length = cigar_op
            if op == 1:
                diff += length
            elif op == 2:
                diff -= length
        tags = []
        tags.append(("XD", int(diff)))
        tags.append(("XS", int(start)))
        tags.append(("XE", int(end)))
        tags.append(("XR", repeat.strip()))
        tags.append(("XC", float(ref)))
        aligned_read.tags = tags
        output_bam.write(aligned_read)
    line = strs.readline()

pysam.sort("temp",output_bam_file)
os.system("mv %s.bam %s"%(output_bam_file, output_bam_file))
pysam.index(output_bam_file)
strs.close()
input_bam.close()
output_bam.close()
os.system("rm temp")

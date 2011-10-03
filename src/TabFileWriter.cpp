/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <cstdlib>
#include <iostream>
#include <list>
#include "TabFileWriter.h"

using namespace std;

TabFileWriter::TabFileWriter(const string& filename):
  TextFileWriter(filename) {
  // write header
  output_stream << "ID\t" <<
    "Nucleotides\t" <<
    "Quality\t"<<
    "MS_start_pos\t" <<
    "MS_end_pos\t" <<
    "MS_repeat_period\t" <<
    "left_flank\t" <<
    "ms_region\t" <<
    "right_flank\t" <<
    "msChrom\t" <<
    "msStart\t"<<
    "msEnd\t"<<
    "msRepeat\t"<<
    "refCopyNum\t"<<
    "lStart\t"<<
    "lEnd\t"<<
    "rStart\t"<<
    "rEnd\t"<<
    "diffFromRef\t"<<
    "reverse\t"<<
    "mismatch\t"<<
    "name"<<
    endl;
}

TabFileWriter::~TabFileWriter(){}

void TabFileWriter::WriteRecord(const MSReadRecord& read) {
  output_stream << read.ID << "\t"
		<< read.nucleotides << "\t"
		<< read.quality_scores << "\t"
		<< read.ms_start << "\t"
		<< read.ms_end << "\t"
		<< read.ms_repeat_best_period << "\t"
		<< read.left_flank_nuc << "\t"
		<< read.detected_ms_region_nuc << "\t"
		<< read.right_flank_nuc << "\t"
		<< read.chrom << "\t"
		<< read.msStart << "\t"
		<< read.msEnd << "\t"
		<< read.msRepeat << "\t"
		<< read.refCopyNum << "\t"
		<< read.lStart << "\t"
		<< read.lEnd << "\t"
		<< read.rStart << "\t"
		<< read.rEnd << "\t"
		<< read.diffFromRef << "\t"
		<< read.reverse << "\t"
		<< read.lDist << "," << read.rDist << "\t"
		<< read.name
		<< endl;
}

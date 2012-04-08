/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <cstdlib>
#include <iostream>
#include <list>
#include "TabFileWriter.h"
#include "runtime_parameters.h"
using namespace std;

TabFileWriter::TabFileWriter(const string& filename):
  TextFileWriter(filename) {
  // write command line argument
  output_stream << user_defined_arguments << endl;
  // write header
  output_stream << "ID" << "\t"
		<< "Nucleotides" << "\t"
		<< "Qualities" << "\t"
		<< "STR region" << "\t"
		<< "Period" << "\t"
		<< "Motif" << "\t"
		<< "chrom" << "\t"
		<< "str_start" << "\t"
		<< "str_end" << "\t"
		<< "read_start" << "\t"
		<< "diffFromRef" << "\t"
		<< "refCopyNum" << "\t"
		<< "reverse" << "\t"
		<< "cigar" << "\t"
		<< "partial" << "\t"
		<< "score";
  if (include_orig_read_start) {
    output_stream << "\tread_start" << endl;
  } else {
    output_stream << endl;
  }
}

TabFileWriter::~TabFileWriter(){}

void TabFileWriter::WriteRecord(const MSReadRecord& read) {
  output_stream << read.ID << "\t"
		<< read.nucleotides << "\t"
		<< read.quality_scores << "\t"
		<< read.detected_ms_region_nuc << "\t"
		<< read.ms_repeat_best_period << "\t"
		<< read.msRepeat << "\t"
		<< read.chrom << "\t"
		<< read.msStart << "\t"
		<< read.msEnd << "\t"
		<< read.read_start << "\t"
		<< read.diffFromRef << "\t"
		<< read.refCopyNum << "\t"
		<< read.reverse << "\t"
		<< read.cigar_string << "\t"
		<< read.partial << "\t"
		<< read.sw_score;
  if (include_orig_read_start) {
    output_stream << "\t" << read.orig_start << endl;
  } else {
    output_stream << endl;
  }
}

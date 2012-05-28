/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include "runtime_parameters.h"
#include "TabFileWriter.h"

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
		<< "score"
		<< "\tmate dist\tstitched" << endl;
}

TabFileWriter::~TabFileWriter(){}

void TabFileWriter::WriteRecord(const ReadPair& read_pair) {
  const int& aligned_read_num = read_pair.aligned_read_num;
  const int& paired_dist = read_pair.treat_as_paired ?
    abs(read_pair.reads.at(aligned_read_num).read_start-
	read_pair.reads.at(1-aligned_read_num).read_start) : -1;
  output_stream << read_pair.reads.at(aligned_read_num).ID << "\t"
		<< read_pair.reads.at(aligned_read_num).nucleotides << "\t"
		<< read_pair.reads.at(aligned_read_num).quality_scores << "\t"
		<< read_pair.reads.at(aligned_read_num).detected_ms_nuc << "\t"
		<< read_pair.reads.at(aligned_read_num).ms_repeat_best_period << "\t"
		<< read_pair.reads.at(aligned_read_num).msRepeat << "\t"
		<< read_pair.reads.at(aligned_read_num).chrom << "\t"
		<< read_pair.reads.at(aligned_read_num).msStart << "\t"
		<< read_pair.reads.at(aligned_read_num).msEnd << "\t"
		<< read_pair.reads.at(aligned_read_num).read_start << "\t"
		<< read_pair.reads.at(aligned_read_num).diffFromRef << "\t"
		<< read_pair.reads.at(aligned_read_num).refCopyNum << "\t"
		<< read_pair.reads.at(aligned_read_num).reverse << "\t"
		<< read_pair.reads.at(aligned_read_num).cigar_string << "\t"
		<< read_pair.reads.at(aligned_read_num).partial << "\t"
		<< read_pair.reads.at(aligned_read_num).sw_score << "\t"
		<< paired_dist << "\t"
		<< (!read_pair.treat_as_paired && paired) << endl;

  if (read_pair.treat_as_paired) {
    output_stream << read_pair.reads.at(1-aligned_read_num).ID << "\t"
		  << read_pair.reads.at(1-aligned_read_num).nucleotides << "\t"
		  << read_pair.reads.at(1-aligned_read_num).quality_scores << "\t"
		  << "-" << "\t"
		  << read_pair.reads.at(aligned_read_num).ms_repeat_best_period << "\t"
		  << read_pair.reads.at(aligned_read_num).msRepeat << "\t"
		  << read_pair.reads.at(aligned_read_num).chrom << "\t"
		  << read_pair.reads.at(aligned_read_num).msStart << "\t"
		  << read_pair.reads.at(aligned_read_num).msEnd << "\t"
		  << read_pair.reads.at(1-aligned_read_num).read_start << "\t"
		  << "-" << "\t"
		  << read_pair.reads.at(aligned_read_num).refCopyNum << "\t"
		  << !read_pair.reads.at(aligned_read_num).reverse << "\t"
		  << read_pair.reads.at(1-aligned_read_num).cigar_string << "\t" 
		  << "-" << "\t"
		  << read_pair.reads.at(1-aligned_read_num).sw_score << "\t"
		  << paired_dist << "\t-1" << endl;
  }
}

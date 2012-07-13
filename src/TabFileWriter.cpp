/*
Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>

This file is part of lobSTR.

lobSTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

lobSTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with lobSTR.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <string>

#include "src/runtime_parameters.h"
#include "src/TabFileWriter.h"

using namespace std;

TabFileWriter::TabFileWriter(const string& filename)
  : TextFileWriter(filename) {
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

TabFileWriter::~TabFileWriter() {}

void TabFileWriter::WriteRecord(const ReadPair& read_pair) {
  const int& aligned_read_num = read_pair.aligned_read_num;
  const int& paired_dist = read_pair.treat_as_paired ?
    abs(read_pair.reads.at(aligned_read_num).read_start-
        read_pair.reads.at(1-aligned_read_num).read_start) : -1;
  output_stream << read_pair.reads.at(aligned_read_num).ID << "\t"
                << read_pair.reads.at(aligned_read_num).nucleotides << "\t"
                << read_pair.reads.at(aligned_read_num).quality_scores << "\t"
                << read_pair.reads.at(aligned_read_num).detected_ms_nuc << "\t"
                << read_pair.reads.at(aligned_read_num).ms_repeat_best_period
                << "\t"
                << read_pair.reads.at(aligned_read_num).repseq << "\t"
                << read_pair.reads.at(aligned_read_num).chrom << "\t"
                << read_pair.reads.at(aligned_read_num).msStart << "\t"
                << read_pair.reads.at(aligned_read_num).msEnd << "\t"
                << read_pair.reads.at(aligned_read_num).read_start << "\t"
                << read_pair.reads.at(aligned_read_num).diffFromRef << "\t"
                << read_pair.reads.at(aligned_read_num).refCopyNum << "\t"
                << read_pair.reads.at(aligned_read_num).reverse << "\t"
                << read_pair.reads.at(aligned_read_num).cigar_string << "\t"
                << read_pair.reads.at(aligned_read_num).partial << "\t"
                << read_pair.reads.at(aligned_read_num).mapq << "\t"
                << paired_dist << "\t"
                << (!read_pair.treat_as_paired && paired) << endl;
  if (read_pair.treat_as_paired) {
    output_stream << read_pair.reads.at(1-aligned_read_num).ID << "\t"
                  << read_pair.reads.at(1-aligned_read_num).nucleotides << "\t"
                  << read_pair.reads.at(1-aligned_read_num).quality_scores
                  << "\t"
                  << "-" << "\t"
                  << read_pair.reads.at(aligned_read_num).ms_repeat_best_period
                  << "\t"
                  << read_pair.reads.at(aligned_read_num).repseq << "\t"
                  << read_pair.reads.at(aligned_read_num).chrom << "\t"
                  << read_pair.reads.at(aligned_read_num).msStart << "\t"
                  << read_pair.reads.at(aligned_read_num).msEnd << "\t"
                  << read_pair.reads.at(1-aligned_read_num).read_start << "\t"
                  << "-" << "\t"
                  << read_pair.reads.at(aligned_read_num).refCopyNum << "\t"
                  << !read_pair.reads.at(aligned_read_num).reverse << "\t"
                  << read_pair.reads.at(1-aligned_read_num).cigar_string
                  << "\t"
                  << "-" << "\t"
                  << read_pair.reads.at(1-aligned_read_num).mapq << "\t"
                  << paired_dist << "\t-1" << endl;
  }
}

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
#include "src/GenotypeTabWriter.h"

using namespace std;

GenotypeTabWriter::GenotypeTabWriter(const string& filename)
  : TextFileWriter(filename) {
  // write command line arguments
  if (user_defined_arguments.size() > 0) {
    output_stream << (user_defined_arguments.
                      substr(4, user_defined_arguments.size()-5)) << endl;
  }
  output_stream << user_defined_arguments_allelotyper << endl;
}

GenotypeTabWriter::~GenotypeTabWriter() {}

void GenotypeTabWriter::WriteRecord(const STRRecord& str_record) {
  if (((str_record.allele1 == -10000 || str_record.allele2 == -10000) &&
        str_record.partial_coverage == 0)
      && (str_record.stop > 0 && str_record.start > 0 &&
          str_record.stop > str_record.start)) {
    return;
  }
  output_stream << str_record.chrom << "\t"
                << str_record.start << "\t"
                << str_record.stop << "\t"
                << str_record.repseq << "\t"
                << str_record.period << "\t"
                << str_record.refcopy << "\t"
                << str_record.allele1_string << ","
                << str_record.allele2_string << "\t"
                << str_record.coverage << "\t"
                << str_record.agreeing << "\t"
                << str_record.coverage - str_record.agreeing << "\t"
                << str_record.readstring << "\t"
                << str_record.score << "\t" // posterior prob of call
                << str_record.allele1_score << "\t" // marginal posterior prob
                << str_record.allele2_score << "\t" // marginal posterior prob
                << str_record.partial_coverage << "\t"
                << str_record.max_partial_string << "\t"
                << str_record.partialreadstring << "\t"
                << str_record.num_stitched << endl;
}

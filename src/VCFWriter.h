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

#ifndef SRC_VCFWRITER_H__
#define SRC_VCFWRITER_H__

#include <stdlib.h>

#include <set>
#include <string>

#include "src/STRRecord.h"
#include "src/TextFileWriter.h"

const int MISSING = -10000;

class VCFWriter : public TextFileWriter {
 public:
  explicit VCFWriter(const std::string& filename,
		     const std::vector<std::string>& samples);
  virtual ~VCFWriter();
  void WriteRecord(const STRRecord& str_record);
  void WriteSample(const STRRecord& str_record, size_t sample_index,
		   std::map<int,int> allele_to_index);
  std::string GetSTRVar(const std::string& refseq,
                        const std::string& ref_repseq,
                        int alleles);
};

#endif  // SRC_VCFWRITER_H__

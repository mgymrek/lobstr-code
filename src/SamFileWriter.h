/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>

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

#ifndef SRC_SAMFILEWRITER_H_
#define SRC_SAMFILEWRITER_H_

#include <string>
#include <map>

#include "api/BamAlignment.h"
#include "api/BamWriter.h"
#include "api/SamHeader.h"
#include "api/SamReadGroup.h"
#include "api/SamReadGroupDictionary.h"
#include "src/ReadPair.h"

namespace BamTools {
  class BamWriter;
  struct BamAlignment;
}

class SamFileWriter {
 public:
  SamFileWriter(const std::string& _filename,
                const std::map<std::string, int>& _chrom_sizes);

  /* Write alignment from lobSTR */
  void WriteRecord(const ReadPair& read_pair);

  /* Write read from allelotype
     - filter referes to filter status. PASS if not filtered
     - chrom, str_start, and str_end refer to specific STRs the read
           was filtered from. Otherwise set to chrom="", str_start=-1, str_end=-1, repseq="", allele=0
   */
  void WriteAllelotypeRead(const BamTools::BamAlignment& aln, const std::string& filter,
                           const std::string& chrom, const int& str_start, const int& str_end,
                           const std::string& repseq, const int& allele);
  virtual ~SamFileWriter();
 private:
  std::string StandardizeReadID(const std::string& readid, bool paired);
  std::map<std::string, int> chrom_sizes;
  BamTools::BamWriter writer;
};

#endif  // SRC_SAMFILEWRITER_H_

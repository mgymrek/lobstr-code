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
  void WriteRecord(const ReadPair& read_pair);
  virtual ~SamFileWriter();
 private:
  std::map<std::string, int> chrom_sizes;
  BamTools::BamWriter writer;
};

#endif  // SRC_SAMFILEWRITER_H_

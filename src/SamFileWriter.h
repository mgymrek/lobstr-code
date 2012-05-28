/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_SAMFILEWRITER_H_
#define SRC_SAMFILEWRITER_H_

#include <string>
#include <map>

#include "api/BamAlignment.h"
#include "api/BamWriter.h"
#include "src/ReadPair.h"

namespace BamTools {
  class BamWriter;
  class BamAlignment;
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

/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_BAMFILEREADER_H__
#define SRC_BAMFILEREADER_H__

#include <string>

#include "src/api/BamReader.h"
#include "src/TextFileReader.h"

namespace BamTools {
  class BamReader;
}

class BamFileReader : public TextFileReader {
 public:
  explicit BamFileReader(const std::string& _filename="");
  virtual bool GetNextRecord(ReadPair* read_pair);
  virtual bool GetNextRead(MSReadRecord* read);
 private:
  BamTools::BamReader reader;
};

#endif  // SRC_BAMFILEREADER_H__

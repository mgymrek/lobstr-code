/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_FASTAFILEREADER_H__
#define SRC_FASTAFILEREADER_H__

#include <string>

#include "src/IFileReader.h"
#include "src/TextFileReader.h"

class FastaFileReader : public TextFileReader {
 public:
  explicit FastaFileReader(const std::string& _filename="");
  virtual bool GetNextRecord(ReadPair* read_pair);
  virtual bool GetNextRead(MSReadRecord* read);
};

#endif  // SRC_FASTAFILEREADER_H__

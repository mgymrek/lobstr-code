/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_FASTQFILEREADER_H__
#define SRC_FASTQFILEREADER_H__

#include <string>

#include "src/TextFileReader.h"

class FastqFileReader : public TextFileReader {
 public:
  explicit FastqFileReader(const std::string& _filename="");
  virtual bool GetNextRecord(ReadPair* read_pair);
  virtual bool GetNextRead(MSReadRecord* read);
};

#endif  // SRC_FASTQFILEREADER_H__

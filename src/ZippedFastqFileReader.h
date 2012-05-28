/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_ZIPPEDFASTQFILEREADER_H__
#define SRC_ZIPPEDFASTQFILEREADER_H__

#include <string>

#include "src/ZippedTextFileReader.h"

class ZippedFastqFileReader : public ZippedTextFileReader {
 public:
  explicit ZippedFastqFileReader(const std::string& _filename="");
  virtual bool GetNextRecord(ReadPair* read_pair);
  virtual bool GetNextRead(MSReadRecord* read);
};

#endif  // SRC_ZIPPEDFASTQFILEREADER_H__

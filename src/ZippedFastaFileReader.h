/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_ZIPPEDFASTAFILEREADER_H__
#define SRC_ZIPPEDFASTAFILEREADER_H__

#include <string>

#include "src/ZippedTextFileReader.h"

class ZippedFastaFileReader : public ZippedTextFileReader {
 public:
  explicit ZippedFastaFileReader(const std::string& _filename="");
  virtual bool GetNextRecord(ReadPair* read_pair);
  virtual bool GetNextRead(MSReadRecord* read);
};

#endif  // SRC_ZIPPEDFASTAFILEREADER_H__

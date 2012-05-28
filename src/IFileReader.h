/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_IFILEREADER_H__
#define SRC_IFILEREADER_H__

#include "src/MSReadRecord.h"
#include "src/ReadPair.h"

class IFileReader {
 public:
  virtual ~IFileReader() { }
  // Get next record as a ReadPair (read2 is NULL in single end mode)
  virtual bool GetNextRecord(ReadPair* read_pair) = 0;
  // Get next read from a file
  virtual bool GetNextRead(MSReadRecord* read) = 0;
};

#endif  // SRC_IFILEREADER_H__

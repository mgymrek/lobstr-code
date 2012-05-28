/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_FASTQPAIREDFILEREADER_H__
#define SRC_FASTQPAIREDFILEREADER_H__

#include <string>

#include "src/FastqFileReader.h"
#include "src/IFileReader.h"

class FastqPairedFileReader : public IFileReader {
 public:
  FastqPairedFileReader(const std::string& _filename1="",
                        const std::string& _filename2="");
  virtual bool GetNextRecord(ReadPair* read_pair);
  virtual bool GetNextRead(MSReadRecord* read);
 private:
  IFileReader* _reader1;
  IFileReader* _reader2;
};

#endif  // SRC_FASTQPAIREDFILEREADER_H__

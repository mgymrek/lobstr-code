/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_FASTAPAIREDFILEREADER_H__
#define SRC_FASTAPAIREDFILEREADER_H__

#include <string>

#include "src/FastaFileReader.h"
#include "src/IFileReader.h"

class FastaPairedFileReader : public IFileReader {
 public:
  FastaPairedFileReader(const std::string& _filename1="",
                        const std::string& _filename2="");
  virtual bool GetNextRecord(ReadPair* read_pair);
  virtual bool GetNextRead(MSReadRecord* read);
 private:
  IFileReader* _reader1;
  IFileReader* _reader2;
};

#endif  // SRC_FASTAPAIREDFILEREADER_H__

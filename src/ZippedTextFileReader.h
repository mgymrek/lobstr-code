/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_ZIPPEDTEXTFILEREADER_H__
#define SRC_ZIPPEDTEXTFILEREADER_H__

#include <string>

#include "src/gzstream.h"
#include "src/IFileReader.h"

class ZippedTextFileReader : public IFileReader {
 public:
  explicit ZippedTextFileReader(const std::string& _filename="");
  virtual ~ZippedTextFileReader();
  virtual bool GetNextRead(MSReadRecord* read);
  virtual bool GetNextRecord(ReadPair* read_pair);
  bool GetNextLine(std::string* line);

 protected:
  size_t current_line;
  std::string filename;
  igzstream *input_file_stream;
  std::istream &input_stream;
  static igzstream* create_file_stream(const std::string& filename);
};

#endif  // SRC_ZIPPEDTEXTFILEREADER_H__

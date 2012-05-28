/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_TEXTFILEWRITER_H_
#define SRC_TEXTFILEWRITER_H_

#include <fstream>
#include <iostream>
#include <string>

#include "src/MSReadRecord.h"

class TextFileWriter {
 public:
  explicit TextFileWriter(const std::string& _filename);
  void WriteRecord(const MSReadRecord& msread);
  void Write(const std::string& str);
  virtual ~TextFileWriter();

 protected:
  std::string filename;
  std::ofstream *output_file_stream;
  std::ostream &output_stream;
  static std::ofstream* create_file_stream(const std::string& filename);
};

#endif  // SRC_TEXTFILEWRITER_H_

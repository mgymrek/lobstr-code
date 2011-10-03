/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __TEXT_FILE_READER_H__
#define __TEXT_FILE_READER_H__

#include <istream>
#include <fstream>
#include <string>

#include "IFileReader.h"

class TextFileReader : public IFileReader {
 public:
  TextFileReader (const std::string& _filename="") ;
  virtual ~TextFileReader();
  virtual bool GetNextRecord(MSReadRecord* read);
  bool GetNextLine(std::string* line);

 protected:
  size_t current_line ;
  std::string filename;
  std::ifstream *input_file_stream;
  std::istream &input_stream;
  static std::ifstream* create_file_stream(const std::string& filename);
};

#endif /* __TEXT_FILE_READER_H__ */

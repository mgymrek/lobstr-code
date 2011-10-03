/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <err.h>
#include <error.h>
#include <iostream>

#include "TextFileWriter.h"

using namespace std;

TextFileWriter::TextFileWriter(const std::string& _filename)
 :filename(_filename),
  output_file_stream((_filename.empty() ? NULL : 
		      create_file_stream(filename))),
  output_stream(*output_file_stream) {}

TextFileWriter::~TextFileWriter() {
  if (output_file_stream!=NULL) {
    delete output_file_stream;
    output_file_stream = NULL ;
  }
}

void TextFileWriter::Write(const string& str) {
  output_stream << str << "\n";
}

std::ofstream* TextFileWriter::create_file_stream(const std::string &filename) {
  std::ofstream *output_stream = new std::ofstream(filename.c_str(), std::ios_base::out);
  if (output_stream==NULL)
    err(1, "Failed to allocate memory for ifstream");
  if ((! (*output_stream)))
    err(1, "Failed to open file '%s'", filename.c_str());
  
  return output_stream;
}

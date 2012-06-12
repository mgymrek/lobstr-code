/*
Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>

This file is part of lobSTR.

lobSTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

lobSTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with lobSTR.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <err.h>

#include <iostream>
#include <string>

#include "src/TextFileWriter.h"

using namespace std;

TextFileWriter::TextFileWriter(const std::string& _filename)
  : filename(_filename),
    output_file_stream((_filename.empty() ? NULL :
                        create_file_stream(filename))),
    output_stream(*output_file_stream) {}

TextFileWriter::~TextFileWriter() {
  if (output_file_stream != NULL) {
    delete output_file_stream;
    output_file_stream = NULL;
  }
}

void TextFileWriter::Write(const string& str) {
  output_stream << str << "\n";
}

std::ofstream* TextFileWriter::create_file_stream(const std::string
                                                  &filename) {
  std::ofstream *output_stream = new std::ofstream(filename.c_str(),
                                                   std::ios_base::out);
  if (output_stream == NULL)
    err(1, "Failed to allocate memory for ifstream");
  if ((!(*output_stream)))
    err(1, "Failed to open file '%s'", filename.c_str());

  return output_stream;
}

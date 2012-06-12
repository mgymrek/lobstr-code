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

#include <fstream>
#include <iostream>
#include <string>

#include "src/IFileReader.h"
#include "src/TextFileReader.h"

using namespace std;

TextFileReader::TextFileReader(const std::string& _filename)
  : current_line(0), filename(_filename),
    input_file_stream(_filename.empty() ? NULL :
                      create_file_stream(filename)),
    input_stream(_filename.empty() ? cin : *input_file_stream ) {}

std::ifstream* TextFileReader::create_file_stream(const std::string
                                                  &filename) {
  std::ifstream *input_stream = new std::ifstream(filename.c_str(),
                                                  std::ios_base::in);
  if (input_stream == NULL)
    err(1, "Failed to allocate memory for ifstream");
  if (!(*input_stream))
    err(1, "Failed to open file '%s'", filename.c_str());
  return input_stream;
}

TextFileReader::~TextFileReader() {
  if (input_file_stream != NULL) {
    delete input_file_stream;
    input_file_stream = NULL;
  }
}

bool TextFileReader::GetNextLine(string* line) {
  current_line++;
  if (!getline(input_stream, *line))
    return false;
  return true;
}

bool TextFileReader::GetNextRecord(ReadPair* read_pair) {
  return false;  // do nothing
}
bool TextFileReader::GetNextRead(MSReadRecord* read) {
  return false;  // do nothing
}

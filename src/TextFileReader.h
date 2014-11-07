/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>

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

#ifndef SRC_TEXTFILEREADER_H__
#define SRC_TEXTFILEREADER_H__

#include <istream>
#include <fstream>
#include <string>

#include "src/IFileReader.h"

class TextFileReader : public IFileReader {
 public:
  explicit TextFileReader(const std::string& _filename="");
  virtual ~TextFileReader();
  virtual bool GetNextRead(MSReadRecord* read);
  virtual bool GetNextRecord(ReadPair* read_pair);
  bool GetNextLine(std::string* line);

 protected:
  size_t current_line;
  std::string filename;
  std::ifstream *input_file_stream;
  std::istream &input_stream;
  static std::ifstream* create_file_stream(const std::string& filename);
};

#endif  // SRC_TEXTFILEREADER_H__

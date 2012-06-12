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

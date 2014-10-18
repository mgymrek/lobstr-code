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

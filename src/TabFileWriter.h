/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_TABFILEWRITER_H__
#define SRC_TABFILEWRITER_H__

#include <stdlib.h>

#include <string>

#include "src/ReadPair.h"
#include "src/TextFileWriter.h"

class TabFileWriter : public TextFileWriter {
 public:
  explicit TabFileWriter(const std::string& filename);
  virtual ~TabFileWriter();
  void WriteRecord(const ReadPair& read_pair);
};

#endif  // SRC_TABFILEWRITER_H__

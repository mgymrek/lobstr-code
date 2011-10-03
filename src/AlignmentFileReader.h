/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef MSTABLEFILEREADER_H_
#define ALIGNMENTTABLEFILEREADER_H_

#include <string>

#include "MSReadRecord.h"
#include "TextFileReader.h"

class AlignmentFileReader : public TextFileReader {
 public:
  AlignmentFileReader(const std::string& _filename);
  virtual ~AlignmentFileReader();
  virtual bool GetNextRecord(MSReadRecord* msread); // do nothing
};

#endif /* ALIGNMENTTABLEFILEREADER_H_ */

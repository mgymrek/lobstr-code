/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __IFILE_READER_H__
#define __IFILE_READER_H__

#include "MSReadRecord.h"

class IFileReader {
 public:
  virtual ~IFileReader() { };
  virtual bool GetNextRecord(MSReadRecord* read) = 0 ;
};

#endif /* __IFILE_READER_H__ */

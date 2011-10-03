/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __FASTQ_FILE_READER_H__
#define __FASTQ_FILE_READER_H__

#include <istream>
#include <fstream>
#include <string>

#include "TextFileReader.h"

class FastqFileReader : public TextFileReader {
 public:
  FastqFileReader ( const std::string& _filename="" ) ;
  virtual bool GetNextRecord(MSReadRecord* read);
};

#endif /*  __FASTQ_FILE_READER_H__ */

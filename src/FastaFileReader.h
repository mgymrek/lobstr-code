/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __FASTA_FILE_READER_H__
#define __FASTA_FILE_READER_H__

#include <istream>
#include <fstream>
#include <string>

#include "IFileReader.h"
#include "TextFileReader.h"

class FastaFileReader : public TextFileReader {
 public:
  FastaFileReader (const std::string& _filename="" ) ;
  virtual bool GetNextRecord(MSReadRecord* read);
};

#endif /* __FASTA_FILE_READER_H__ */

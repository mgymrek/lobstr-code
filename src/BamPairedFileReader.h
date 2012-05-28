/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __BAM_PAIRED_FILE_READER_H__
#define __BAM_PAIRED_FILE_READER_H__

#include <istream>
#include <fstream>
#include <string>

#include "api/BamReader.h"
#include "TextFileReader.h"

namespace BamTools {
  class BamReader;
}

class BamPairedFileReader : public TextFileReader {
 public:
  BamPairedFileReader ( const std::string& _filename="" ) ;
  virtual bool GetNextRecord(ReadPair* read_pair);
  virtual bool GetNextRead(MSReadRecord* read);
 private:
  BamTools::BamReader reader;
};

#endif /*  __BAM_PAIRED_FILE_READER_H__ */

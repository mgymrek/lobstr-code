/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SAMFILEWRITER_H_
#define SAMFILEWRITER_H_

#include <string>
#include <map>

#include "api/BamAlignment.h"
#include "api/BamWriter.h"
#include "MSReadRecord.h"

using namespace std;

namespace BamTools {
  class BamWriter;
  class BamAlignment;
}

class SamFileWriter {
 public:
  SamFileWriter(const std::string& _filename,
		const std::map<string, int>& _chrom_sizes);
  void WriteRecord(const MSReadRecord& msread);
  void WriteAdjustedRecord(const MSReadRecord& msread);
  virtual ~SamFileWriter();
 private:
  std::map<string, int> chrom_sizes;
  BamTools::BamWriter writer;
  string qname;
  int flag;
  string rname;
  int lpos;
  int rpos;
  int mapq;
  string mrnm;
  int mpos;
  int isize;
  string lseq;
  string middle;
  string rseq;
  string cigar;
  string read;
  string qual;
};

#endif /* SAMFILEWRITER_H_ */

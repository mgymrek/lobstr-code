/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef MSTABLEFILEREADER_H_
#define MSTABLEFILEREADER_H_

#include <string>

#include "GenomeReader.h"
#include "MSReadRecord.h"
#include "MSRecord.h"
#include "TextFileReader.h"

class MSTableFileReader : public TextFileReader {
 public:
  MSTableFileReader(const std::string& _filename,
		    const std::string& _genome_filename);
  bool GetNextRecord(MSRecord* msrec);
  virtual ~MSTableFileReader();
  virtual bool GetNextRecord(MSReadRecord* msread); // do nothing
  // get the size of each chromosome
  void GetChromSizes(std::map<string, int>* chrom_sizes);
 private:
  GenomeReader* gReader;
};

#endif /* MSTABLEFILEREADER_H_ */

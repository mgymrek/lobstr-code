/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef GENOMEREADER_H_
#define GENOMEREADER_H_

#include <istream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "TextFileReader.h"

using namespace std;

class GenomeReader: public TextFileReader {
 public:
  GenomeReader(const std::string& _filename);
  virtual ~GenomeReader();
  bool GetGenomeCoords(const std::string& chrom, int start,
		       int end, std::string* nucs);

  // does nothing, placeholder for virtual function so we can extend
  // TextFileReader
  virtual bool GetNextRecord(MSReadRecord* read);

  // map of chr->nucleotides
  map<std::string, std::string> genome;

  // map of chr-> length
  map<std::string, int> lengths;
};

#endif /* GENOMEREADER_H_ */

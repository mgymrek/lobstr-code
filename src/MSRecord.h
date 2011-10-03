/*
  Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef MSRECORD_H_
#define MSRECORD_H_

#include <string>

class MSRecord {
 public:
  // id of this STR record
  int seqid;
  // chromosome of the STR record
  std::string chrom;
  // start location of the STR
  int start;
  // end location of the STR
  int end;
  // how many nucleotides the reference was extended
  int extend;
  // repeated STR sequence
  std::string repeat;
  // copy number in reference
  float copynum;
  // left flank in reference
  std::string leftFlank;
  // right flank in reference
  std::string rightFlank;
  // name of this locus
  std::string name;
};

#endif /* MSRECORD_H_ */

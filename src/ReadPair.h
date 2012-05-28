/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_READPAIR_H_
#define SRC_READPAIR_H_

#include <vector>

#include "src/MSReadRecord.h"

/* 
   A ReadPair object contains information about
   one set of paired end reads. If lobSTR is not running
   in paired end mode, this object is still used, and
   the second read of the pair is set to null.
*/

class ReadPair {
 public:
  // Read count
  size_t read_count;
  // Reads in a pair
  std::vector<MSReadRecord> reads;

  // First read passed detection
  bool read1_passed_detection;
  // Second read passed detection
  bool read2_passed_detection;
  // First read passed alignment
  bool read1_passed_alignment;
  // Second read passed alignment
  bool read2_passed_alignment;
  // Found a unique good alignment
  bool found_unique_alignment;
  // Number of the aligned read
  int aligned_read_num;
  // Treat the output as paired end
  bool treat_as_paired;
};

#endif  // SRC_READPAIR_H_

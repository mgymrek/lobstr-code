/*
Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>

This file is part of lobSTR.

lobSTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

lobSTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with lobSTR.  If not, see <http://www.gnu.org/licenses/>.
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
  // String for alternate mappings
  std::string alternate_mappings;
};

#endif  // SRC_READPAIR_H_

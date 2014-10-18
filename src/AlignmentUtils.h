/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>

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

#ifndef SRC_ALIGNMENTUTILS_H_
#define SRC_ALIGNMENTUTILS_H_

#include <string>

#include "src/Alignment.h"
#include "src/cigar.h"
#include "src/MSReadRecord.h"
#include "src/ReadPair.h"

namespace AlignmentUtils {

  // Calculate map quality score
  int GetMapq(const std::string& aligned_sw_string,
              const std::string& ref_sw_string,
              const std::string& aligned_quals,
              int* edit_dist);

  // Try stitching a pair of reads.
  // Update info in num_aligned_read and
  // treat as single alignment
  bool StitchReads(ReadPair* read_pair,
                   ALIGNMENT* left_alignment,
                   ALIGNMENT* right_alignment);

  // Refine the cigar score and recalculate number of repeats
  bool GetSTRAllele(MSReadRecord* aligned_read,
                    const CIGAR_LIST& cigar_list);
}

#endif  // SRC_ALIGNMENTUTILS_H_

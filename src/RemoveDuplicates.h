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

#ifndef SRC_REMOVEDUPLICATES_H_
#define SRC_REMOVEDUPLICATES_H_

#include <list>

#include "src/AlignedRead.h"

namespace RemoveDuplicates {
  /*
    Remove PCR duplicates
    Reads are duplicates if they have the same start coordinate and same length.
    Choose as representative read the read with the highest quality score
  */
  void RemovePCRDuplicates(std::list<AlignedRead>* aligned_reads);

  /* Get representative read from group of duplicates, read with highest qual score */
  void GetRepRead(const std::list<AlignedRead>& aligned_reads, AlignedRead* rep_alignment);

  /* Get quality score of a read */
  float GetScore(const std::string& quality_string);
}

#endif  // SRC_REMOVEDUPLICATES_H_

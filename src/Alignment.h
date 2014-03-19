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

#ifndef SRC_ALIGNMENT_H_
#define SRC_ALIGNMENT_H_

// Store a single alignment found by BWA
struct ALIGNMENT {
  // Identifier of the STR aligned to
  int id;
  // Chrom of the STR aligned to
  std::string chrom;
  // Start of the STR
  int start;
  // End of the STR
  int end;
  // This alignment is on the left side
  bool left;
  // Repeat motif
  std::string repeat;
  // true = positive, false = minus
  bool strand;
  // Position of the start of the alignment
  int pos;
  // Position of the end of the alignment
  int endpos;
  // Copy number
  float copynum;
  // provide overloaded operators to compare
  // when used as a map key
  bool operator<(const ALIGNMENT& ref_pos1) const {
    if (chrom != ref_pos1.chrom) return true;
    if (start == ref_pos1.start) {
      return end < ref_pos1.end;
    } else {
      return start < ref_pos1.start;
    }
  }
};

#endif  // SRC_ALIGNMENT_H_

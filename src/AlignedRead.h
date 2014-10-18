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

#ifndef SRC_ALIGNEDREAD_H_
#define SRC_ALIGNEDREAD_H_

#include "src/api/BamAlignment.h"

/*
  Struct to keep track of aligned read used in allelotyping
 */
struct AlignedRead {
  std::string ID;
  std::string chrom;
  int msStart;
  int msEnd;
  std::string read_group; // identifies a unique sample
  int read_start;
  std::string nucleotides;
  std::string qualities;
  std::vector<BamTools::CigarOp> cigar_ops;
  std::string repseq;
  int period;
  int diffFromRef;
  float refCopyNum;
  int mate;
  bool strand;
  int stitched;
  int matedist;
  int mapq;
  bool stutter;
};
#endif  // SRC_ALIGNEDREAD_H_

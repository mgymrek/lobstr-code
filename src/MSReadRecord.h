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

#ifndef SRC_MSREADRECORD_H_
#define SRC_MSREADRECORD_H_

#include <string>
#include <vector>

#include "src/cigar.h"

/*
  MSReadRecord stores all information about a single read.
  This may be a single read, or one read of a pair.
  Paired reads are containd in a ReadPair object.
 */

class MSReadRecord {
 public:
  // *** input data *** //
  // ID of the read
  std::string ID;
  // nucleotides of the read
  std::string nucleotides;
  // quality scores of the read
  std::string quality_scores;
  // keep original nucleotides and qualities for printing at the end
  std::string orig_nucleotides;
  std::string orig_qual;

  // *** set during detection *** //
  int ms_start;
  // location in the read whre the STR ends
  int ms_end;
  // left flanking nucleotides
  std::string left_flank_nuc;
  int left_flank_index_from_start;
  // STR nucleotides (detected by STRDetector.cpp)
  std::string detected_ms_region_nuc;
  // right flanking nucleotides
  std::string  right_flank_nuc;
  int right_flank_index_from_end;

  // *** set during alignment *** //
  // Repeat motif
  std::string repseq;
  // start of the original read (before trimming)
  // Note: this field is deprecated and was used only
  // for coordination with STRangers
  int orig_start;
  // reverse complement of ther ead was aligned
  bool reverse;
  // chromosome of STR
  std::string chrom;
  // ID of str
  int strid;
  // start position of STR
  int msStart;
  // end position of STR
  int msEnd;
  // reference copy number
  float refCopyNum;
  // start coord of left flank
  int lStart;
  // end coord of left flank
  int lEnd;
  // start coord of right flank
  int rStart;
  // start coord of end flank
  int rEnd;
  // difference in length from reference
  int diffFromRef;
  // name of this STR locus
  std::string name;
  // alignment score
  int mapq;
  // edit distance from reference, for SAM NM tag
  int edit_dist;
  // cigar score for sam format
  std::vector<CIGAR> cigar;
  // cigar score as a string
  std::string cigar_string;
  // three prime pos
  int read_start;
  // five prime pos
  int read_end;
  // detected STR (reset by BWAReadAligner.cpp)
  std::string  detected_ms_nuc;
  // is the read paired?
  bool paired;
};

#endif  // SRC_MSREADRECORD_H_

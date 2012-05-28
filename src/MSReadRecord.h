/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
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
  // detected rep seq
  std::string repseq;
  // is canonical repseq the reverse of how it appears?
  bool repseq_reverse;
  // location in the read where the STR starts
  int ms_start;
  // location in the read whre the STR ends
  int ms_end;
  // best period determined for the STR
  int ms_repeat_best_period;
  // next best period determined for the STR
  int ms_repeat_next_best_period;
  // left flanking nucleotides
  std::string left_flank_nuc;
  int left_flank_index_from_start;
  // STR nucleotides (detected by STRDetector.cpp)
  std::string detected_ms_region_nuc;
  // right flanking nucleotides
  std::string  right_flank_nuc;
  int right_flank_index_from_end;

  // *** set during alignment *** //
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
  // repeat of aligned read
  std::string msRepeat;
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
  // partially covered
  bool partial;
  // name of this STR locus
  std::string name;
  // alignment score
  int sw_score;
  // cigar score for sam format
  std::vector<CIGAR> cigar;
  // cigar score as a string
  std::string cigar_string;
  // three prime pos
  int read_start;
  // five prime pos
  int read_end;
  // is the left flank all repeats
  bool left_all_repeats;
  // is the right flank all repeats
  bool right_all_repeats;
  // detected STR (reset by BWAReadAligner.cpp)
  std::string  detected_ms_nuc;
};

#endif  // SRC_MSREADRECORD_H_

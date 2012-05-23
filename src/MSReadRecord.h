/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef MSREADRECORD_H_
#define MSREADRECORD_H_

#include <iostream>
#include <istream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "cigar.h"

class MSReadRecord {
public:
  size_t read_counter;
  // input data
  // ID of the read
  std::string ID;
  // nucleotides of the read
  std::string nucleotides;
  // quality scores of the read
  std::string quality_scores;
  // keep original nucleotides and qualities for printing at the end
  std::string orig_nucleotides;
  std::string orig_qual;
  int orig_start;

  // detected rep seq
  std::string repseq;
  // is canonical repseq the reverse of how it appears?
  bool repseq_reverse;

  // detection data
  // location in the read where the STR starts
  int	ms_start;
  // location in the read whre the STR ends
  int	ms_end;
  // best period determined for the STR
  int	ms_repeat_best_period;
  // next best period determined for the STR
  int ms_repeat_next_best_period;
  // left flanking nucleotides
  std::string	left_flank_nuc;
  int left_flank_index_from_start;
  // STR nucleotides (detected by STRDetector.cpp)
  std::string	detected_ms_region_nuc;
  // right flanking nucleotides
  std::string  right_flank_nuc;
  int right_flank_index_from_end;
  // detected STR (reset by BWAReadAligner.cpp)
  std::string  detected_ms_nuc;
  
  // reverse complement of ther ead was aligned
  bool reverse;
  
  // output data
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
  // mismatches in left flank
  int lDist;
  // micmatches in right flank
  int rDist;
  // name of this STR locus
  std::string name;
  // alignment score
  int sw_score;
  // cigar score for sam format
  std::vector<CIGAR> cigar;
  std::string cigar_string;
  // three prime pos
  int read_start;
  // five prime pos
  int read_end;
  // did we say this read is aligned?
  bool accepted;
};

#endif /* MSREADRECORD_H_ */

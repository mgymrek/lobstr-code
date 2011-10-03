/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef READRECORD_H_
#define READRECORD_H_

#include <string>

class ReadRecord {
public:
  size_t read_counter;
  // input data
  // ID of the read
  std::string ID;
  // nucleotides of the read
  std::string nucleotides;
  // quality scores of the read
  std::string quality_scores;

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
  // STR nucleotides
  std::string	detected_ms_region_nuc;
  // right flanking nucleotides
  std::string  right_flank_nuc;
  // detected STR
  std::string  detected_ms_nuc;
  
  // reverse complement of ther ead was aligned
  bool reverse;
  
  // output data
  // chromosome of aligned read
  std::string chrom;
  // start position of aligned read
  int msStart;
  // end position of aligned read
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
  // mismatches in left flank
  int lDist;
  // micmatches in right flank
  int rDist;
  // name of this STR locus
  std::string name;
};

#endif /* MSREADRECORD_H_ */

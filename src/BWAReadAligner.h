/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef BWA_READALIGNER_H_
#define BWA_READALIGNER_H_

#include <map>
#include <set>
#include <vector>
#include <list>

#include "bwtaln.h"
#include "bwase.h"
#include "common.h"
#include "MSReadRecord.h"
#include "MSRecord.h"
#include "SamFileWriter.h"
#include "TabFileWriter.h"


struct ALIGNMENT {
  int id;
  bool left;
  std::string chrom;
  int start;
  int end;
  std::string repeat;
  float copynum;
  std::string name;
  bool strand;
  int pos;
  int score;
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

class BWAReadAligner {
 public:
  BWAReadAligner(map<std::string, BWT>* bwt_references,
		 map<std::string, BNT>* bnt_annotations,
		 gap_opt_t *opts);
  virtual ~BWAReadAligner();

  // main function - call BWA aligner 
  bool ProcessRead(MSReadRecord* read);

 protected:
  void ParseRefid(const std::string& refstring, ALIGNMENT* refid);
  
  bool GetAlignmentCoordinates(bwa_seq_t* aligned_seqs, const std::string& repseq,
			       list<ALIGNMENT>* alignments);

  static bool GetSharedAln(const list<ALIGNMENT>& map1,
			   const list<ALIGNMENT>& map2,
			   ALIGNMENT* left_refid,
			   ALIGNMENT* right_refid);

  map<std::string, BWT>* _bwt_references;
  map<std::string, BNT>* _bnt_annotations;
  gap_opt_t *_opts;
};

#endif /* BWA_READALIGNER_H_ */

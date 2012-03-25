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
#include "TabFileWriter.h"

using namespace std;

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
		 map<int, REFSEQ>* ref_sequences,
		 gap_opt_t *opts);
  virtual ~BWAReadAligner();

  // main function - call BWA aligner 
  bool ProcessRead(MSReadRecord* read);

 protected:
  void ParseRefid(const std::string& refstring, ALIGNMENT* refid);
  
  bool GetAlignmentCoordinates(bwa_seq_t* aligned_seqs, const std::string& repseq,
			       std::list<ALIGNMENT>* alignments);

  static bool GetSharedAln(const list<ALIGNMENT>& map1,
			   const list<ALIGNMENT>& map2,
			   ALIGNMENT* left_refid,
			   ALIGNMENT* right_refid);

  // Perform local realignment, adjust exact STR boundaries
  // update cigar score.
  bool AdjustAlignment(MSReadRecord* aligned_read, bool partial, bool left_aligned, bool right_aligned);

  // Get the number of repeat units using the adjusted cigar score
  // Also refine position of flanking regions
  bool GetSTRAllele(MSReadRecord* aligned_read, const CIGAR_LIST& cigar_list, bool* partial);

  map<std::string, BWT>* _bwt_references;
  map<std::string, BNT>* _bnt_annotations;
  map<int, REFSEQ>* _ref_sequences;
  gap_opt_t *_opts;
};

#endif /* BWA_READALIGNER_H_ */

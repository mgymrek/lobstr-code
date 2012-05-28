/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdexcept>

#include "BWAReadAligner.h"
#include "bwaseqio.h"
#include "nw.h"
#include "runtime_parameters.h"

using namespace std;

extern unsigned char nst_nt4_table[256];
// Number of N's used to pad each reference
int PAD = 50;
// max size of cigar score to allow
// more than this is likely bad alignment
size_t MAX_CIGAR_SIZE = 6;
// For partial alignment, clear flanking region
// alignments if maps all over, likely
// repetitive
size_t MAX_ALLOWED_FLANK = 10;
// Maximum difference between mate alignment
// and STR read alignment
int MAX_PAIRED_DIFF = 1000;
// BWA % mismatching for mate alignment
float MATE_FNR = 0.06;
// Minimum number of bp for stitch overlap
size_t MIN_STITCH_OVERLAP = 16;
// Percent identity required to stitch
float STITCH_REQUIRED_SCORE = 0.8;
// Allowed difference in score between returned stitch
// and next best stitch
float STITCH_DIFF = 0;

BWAReadAligner::BWAReadAligner(map<std::string, BWT>* bwt_references,
			       map<std::string, BNT>* bnt_annotations,
			       map<int, REFSEQ>* ref_sequences,
			       gap_opt_t *opts) {
  bwase_initialize();
  _bwt_references = bwt_references;
  _bnt_annotations = bnt_annotations;
  _ref_sequences = ref_sequences;
  _opts = opts;
  _default_opts = gap_init_opt();
  _default_opts->fnr = MATE_FNR;
}

bool BWAReadAligner::ProcessReadPair(ReadPair* read_pair) {
  // Initialize status variables
  read_pair->read1_passed_alignment = false;
  read_pair->read2_passed_alignment = false;
  read_pair->found_unique_alignment = false;
  read_pair->aligned_read_num = -1;

  if (paired) {
    // all valid alignments for individual reads in the pair
    vector<ALIGNMENT> good_left_alignments_read1;
    vector<ALIGNMENT> good_right_alignments_read1;
    vector<ALIGNMENT> good_left_alignments_read2;
    vector<ALIGNMENT> good_right_alignments_read2;
    
    // Step 1: Align each read separately
    if (read_pair->read1_passed_detection) {
      if (ProcessRead(&read_pair->reads.at(0),
                      &good_left_alignments_read1,
                      &good_right_alignments_read1)) {
        read_pair->read1_passed_alignment = true;
      }
    }
    if (read_pair->read2_passed_detection) {
      if (ProcessRead(&read_pair->reads.at(1),
                      &good_left_alignments_read2,
                      &good_right_alignments_read2)) {
        read_pair->read2_passed_alignment = true;
      }
    }
    // Need at least one read to pass
    if (!(read_pair->read1_passed_alignment ||
          read_pair->read2_passed_alignment)) {
      return false;
    }
    // Step 2: Determine if unique valid alignment
    // NOTE: for now, assume only one of the pair aligns
    if (read_pair->read1_passed_alignment ||
        read_pair->read2_passed_alignment) {
      ALIGNMENT matealign;
      size_t num_alignments = 0;
      size_t index_of_hit;
      // Get info for the read that aligned
      if (read_pair->read1_passed_alignment) {
        read_pair->aligned_read_num = 0;
        num_alignments = good_left_alignments_read1.size();
      } else {
        read_pair->aligned_read_num = 1;
        num_alignments = good_left_alignments_read2.size();
      }
      const vector<ALIGNMENT>& good_left = 
        read_pair->aligned_read_num == 0 ?
        good_left_alignments_read1 :
        good_left_alignments_read2;
      const vector<ALIGNMENT>& good_right =
        read_pair->aligned_read_num == 0 ?
        good_right_alignments_read1 :
        good_right_alignments_read2;
      
      // For each alignment, check other read
      vector<ALIGNMENT> mate_alignments;
      if (!AlignMate(*read_pair, &mate_alignments,
                     read_pair->reads.at(read_pair->aligned_read_num).repseq)) {
        return false;
      }
      for (size_t i = 0; i < num_alignments; i++) {
        const ALIGNMENT& lalign = good_left.at(i);
        const ALIGNMENT& ralign = good_right.at(i);
        if (CheckMateAlignment(mate_alignments, lalign, ralign,
                               &matealign)) {
          if (!read_pair->found_unique_alignment) {
            index_of_hit = i;
            read_pair->found_unique_alignment = true;
          } else {
            return false; // multiple mapper
          }
        }
      }
      
      // Step 3: Adjust alignment and output
      if (read_pair->found_unique_alignment) {
        ALIGNMENT final_left_alignment =
          good_left.at(index_of_hit);
        ALIGNMENT final_right_alignment =
          good_right.at(index_of_hit);
        // try stitching first
        if (StitchReads(read_pair, &final_left_alignment,
                        &final_right_alignment)) {
          return OutputAlignment(read_pair, final_left_alignment,
                                 final_right_alignment,
                                 matealign, false);
        } else {
          return OutputAlignment(read_pair, final_left_alignment,
                                 final_right_alignment,
                                 matealign, true);
        }
      }
    }
  } else { // single end reads
    // all valid alignments for each read
    vector<ALIGNMENT> good_left_alignments_read1;
    vector<ALIGNMENT> good_right_alignments_read1;
    if (read_pair->read1_passed_detection) {
      if (ProcessRead(&read_pair->reads.at(0),
                      &good_left_alignments_read1,
                      &good_right_alignments_read1)) {
        
        // Get rid of multi mappers
        if (good_left_alignments_read1.size() > 1) {
          return false;
        }
        
        // Output the alignment
        read_pair->read1_passed_alignment = true;
        read_pair->found_unique_alignment = true;
        ALIGNMENT matealign; // dummy
        return OutputAlignment(read_pair, good_left_alignments_read1.front(),
                               good_right_alignments_read1.front(),
                               matealign, false);
      }
    }
  }
  return false;
}

bool BWAReadAligner::ProcessRead(MSReadRecord* read,
				 vector<ALIGNMENT>* good_left_alignments,
				 vector<ALIGNMENT>* good_right_alignments) {
  // do the flanks consist of all repeats?
  bool left_all_repeats = false;
  bool right_all_repeats = false;
  read->partial = false;

  // Detected motif
  const string& repseq = read->repseq;

  // check that we have that repseq
  if (_bwt_references->find(repseq) == _bwt_references->end()) {
    return false;
  }
  
  // Check if flanks are perfect repeats
  if (IsPerfectRepeat(read->left_flank_nuc, repseq)) 
    left_all_repeats = true;
  if (IsPerfectRepeat(read->right_flank_nuc, repseq))
    right_all_repeats = true;
  // don't continue if both are fully repetitive
  if (left_all_repeats && right_all_repeats) return false;

  // Align the flanking regions
  bwa_seq_t* seqs = BWAAlignFlanks(*read);
  bwa_seq_t* seq_left = &seqs[0];
  bwa_seq_t* seq_right = &seqs[1];

  // fill in alignment coordinates
  vector<ALIGNMENT> left_alignments;
  vector<ALIGNMENT> right_alignments;
  if (!left_all_repeats) {
    if (!GetAlignmentCoordinates(seq_left, repseq, &left_alignments)) {
      bwa_free_read_seq(2, seqs);
      return false;
    }
  }
  if (!right_all_repeats) {
    if (!GetAlignmentCoordinates(seq_right, repseq, &right_alignments)) {
      bwa_free_read_seq(2, seqs);
      return false;
    }
  }

  // don't need seqs anymore
  bwa_free_read_seq(2, seqs);

  // get shared alignments
  // If none found, check for partial coverage
  vector<ALIGNMENT> left_refids;
  vector<ALIGNMENT> right_refids;
  if (!GetSharedAlns(left_alignments, right_alignments,
		     &left_refids, &right_refids)) {
    left_refids.clear();
    right_refids.clear();
    // if one flank maps everywhere, clear it
    bool left_cleared = false;
    bool right_cleared = false;
    if (left_alignments.size() > MAX_ALLOWED_FLANK) {
      left_alignments.clear();
      left_cleared = true;
    }
    if (right_alignments.size() > MAX_ALLOWED_FLANK) {
      right_alignments.clear();
      right_cleared = true;
    }
    // check if one flanking region aligns uniquely
    if (left_alignments.size() + right_alignments.size() == 1) {
      read->partial = true;
      // Case 1: anchored the left side
      if (left_alignments.size() != 0 && 
	  (right_all_repeats||right_cleared)) {
	right_all_repeats = true;
	// left side is anchored
	const ALIGNMENT& left_refid =
	  left_alignments.front();
	// update right side
	ALIGNMENT right_refid;
	UpdatePartialFlank(left_refid, &right_refid,
			   (int)read->nucleotides.length(),
			   (int)read->left_flank_nuc.length(),
			   (int)read->right_flank_nuc.length());
	left_refids.push_back(left_refid);
	right_refids.push_back(right_refid);
      }
      // Case 2: anchored the right side
      else if (right_alignments.size() != 0 && 
	       (left_all_repeats||left_cleared)) {
	left_all_repeats = true;
	// Right side is anchored
	const ALIGNMENT& right_refid =
	  right_alignments.front();
	// update fields in left
	ALIGNMENT left_refid;
	UpdatePartialFlank(right_refid, &left_refid,
			   (int)read->nucleotides.length(),
			   (int)read->right_flank_nuc.length(),
			   (int)read->left_flank_nuc.length());
	left_refids.push_back(left_refid);
	right_refids.push_back(right_refid);
      } else {
	return false;
      }
    } else {
      return false;
    }
  }

  read->left_all_repeats = left_all_repeats;
  read->right_all_repeats = right_all_repeats;
  *good_left_alignments = left_refids;
  *good_right_alignments = right_refids;

  // Return true if at least one good alignment
  return (good_left_alignments->size() >= 1 &&
	  good_right_alignments->size() >= 1);
}

bwa_seq_t* BWAReadAligner::BWAAlignFlanks(const MSReadRecord& read) {
  // set up bwa
  int is_comp = _opts->mode&BWA_MODE_COMPREAD;
  bwa_seq_t *seqs, *seq_left, *seq_right;
  seqs=(bwa_seq_t*)calloc(2, sizeof(bwa_seq_t));
  seq_left = &seqs[0];
  seq_right= &seqs[1];

  // left in for historical reasons since BWA reverses sequence
  const string& left_flank_nuc = reverseComplement(read.left_flank_nuc);
  const string& right_flank_nuc = reverseComplement(read.right_flank_nuc);

  // LEFT FLANKING REGION
  const string& left_qual =  reverse(read.quality_scores.
				     substr(0, left_flank_nuc.length()));
  seq_left->bc[0] = 0;
  seq_left->tid = -1;
  seq_left->qual = 0;
  seq_left->full_len = seq_left->clip_len = seq_left->len =
    left_flank_nuc.length();
  seq_left->seq = (ubyte_t*)calloc(seq_left->len, 1);
  seq_left->qual = (ubyte_t*)calloc(seq_left->len, 1);
  for (int i = 0; i != seq_left->full_len; ++i) {
    seq_left->seq[i] = nst_nt4_table[(int)left_flank_nuc.at(i)];
    if (fastq) {
      seq_left->qual[i] = left_qual.at(i)+33<126?left_qual.at(i)+33:126;
    } else {
      seq_left->qual[i] = 126;
    }
  }
  seq_left->rseq = (ubyte_t*)calloc(seq_left->full_len,1);
  memcpy(seq_left->rseq, seq_left->seq, seq_left->len);
  seq_reverse(seq_left->len, seq_left->seq, 0);
  seq_reverse(seq_left->len, seq_left->rseq, is_comp);
  seq_reverse(seq_left->len, seq_left->qual,0);
  seq_left->name = strdup((const char*)read.ID.c_str());
				   
  // RIGHT FLANKING REGION
  const string& right_qual = reverse(read.quality_scores.
    substr(read.quality_scores.size() - right_flank_nuc.length(),
	   right_flank_nuc.length()));
  seq_right->bc[0] = 0;
  seq_right->tid = -1;
  seq_right->qual = 0;
  seq_right->full_len = seq_right->clip_len = seq_right->len =
    right_flank_nuc.length();
  seq_right->seq = (ubyte_t*)calloc(seq_right->len, 1);
  seq_right->qual = (ubyte_t*)calloc(seq_right->len, 1);
  for (int i = 0; i != seq_right->full_len; ++i) {
    seq_right->seq[i] = nst_nt4_table[(int)right_flank_nuc.at(i)];
    if (fastq) {
      seq_right->qual[i] = right_qual.at(i)+33<126?right_qual.at(i)+33:126;
    } else {
      seq_right->qual[i] = 126;
    }
  }
  seq_right->rseq = (ubyte_t*)calloc(seq_right->full_len,1);
  memcpy(seq_right->rseq, seq_right->seq, seq_right->len);
  seq_reverse(seq_right->len, seq_right->seq, 0);
  seq_reverse(seq_right->len, seq_right->rseq, is_comp);
  seq_reverse(seq_right->len, seq_right->qual, 0);
  seq_right->name = strdup((const char*)read.ID.c_str());

  // call bwa with appropriate options
  bwa_cal_sa_reg_gap(0, _bwt_references->at(read.repseq).bwt, 2, seqs, _opts); 

  return seqs;
}

void BWAReadAligner::UpdatePartialFlank(const ALIGNMENT& aligned_flank,
					ALIGNMENT* unaligned_flank,
					int nuc_size, int aligned_flank_size,
					int unaligned_flank_size) {
  unaligned_flank->left = !aligned_flank.left;
  unaligned_flank->chrom = aligned_flank.chrom;
  unaligned_flank->start = aligned_flank.start - extend;
  unaligned_flank->end = aligned_flank.end - extend;
  unaligned_flank->id = aligned_flank.id;
  unaligned_flank->copynum = aligned_flank.copynum;
  unaligned_flank->repeat = aligned_flank.repeat;
  if (!unaligned_flank->left) {
    unaligned_flank->pos = aligned_flank.pos +
      nuc_size - unaligned_flank_size;
  } else {
    unaligned_flank->pos = aligned_flank.pos +
      + aligned_flank_size - nuc_size;
  }
}

void BWAReadAligner::ParseRefid(const string& refstring, ALIGNMENT* refid) {
  if (_ref_info.find(refstring) != _ref_info.end()) {
    const ALIGNMENT& ref_info = _ref_info.at(refstring);
    refid->id = ref_info.id;
    refid->left = ref_info.left;
    refid->chrom = ref_info.chrom;
    refid->start = ref_info.start;
    refid->end = ref_info.end;
    refid->copynum = ref_info.copynum;
    refid->name = ref_info.name;
  } else { // parse and store
    ALIGNMENT ref_info;
    vector<string> items;
    split(refstring, '_', items);
    refid->id = atoi(items.at(0).c_str());
    refid->left = (items.at(1) == "L");
    refid->chrom = items.at(2);
    refid->start = atoi(items.at(3).c_str())+extend;
    refid->end = atoi(items.at(4).c_str());
    refid->copynum = atof(items.at(6).c_str());
    if (items.size() > 7)
      refid->name = items.at(7);
    // update cache
    ref_info.id = refid->id;
    ref_info.left = refid->left;
    ref_info.chrom = refid->chrom;
    ref_info.start = refid->start;
    ref_info.end = refid->end;
    ref_info.copynum = refid->copynum;
    ref_info.name = refid->name;
    _ref_info.insert(pair<string,ALIGNMENT>(refstring, ref_info));
  }
}

// ** copied from BWA ** //
int64_t pos_end_multi(const bwt_multi1_t *p, int len) // analogy to pos_end()
{
	if (p->cigar) {
		int j;
		int64_t x = p->pos;
		for (j = 0; j != p->n_cigar; ++j) {
			int op = __cigar_op(p->cigar[j]);
			if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
		}
		return x;
	} else return p->pos + len;
}


bool BWAReadAligner::GetAlignmentCoordinates(bwa_seq_t* aligned_seqs,
					     const std::string&repseq,
					     vector<ALIGNMENT>* alignments) {
  // fill in alignment properties
  bwa_aln2seq_core(aligned_seqs->n_aln, aligned_seqs->aln, aligned_seqs,0,1000);
  bwa_cal_pac_pos_core(0, _bwt_references->at(repseq).bwt[1-aligned_seqs->strand],
		       &aligned_seqs[0], _opts->max_diff, _opts->fnr);

  // get coords
  for (int i = 0; i < aligned_seqs->n_multi; ++i) {
    bwt_multi1_t *q = aligned_seqs->multi + i;
    // set pos
    if (!q->strand) {
      q->pos = _bwt_references->at(repseq).bwt[1]->seq_len -
	bwt_sa(_bwt_references->at(repseq).bwt[1], q->pos) +
	aligned_seqs[0].len;
    } else {
      q->pos = bwt_sa(_bwt_references->at(repseq).bwt[0], q->pos);
    }

    int j, seqid; 
    j = pos_end_multi(q, aligned_seqs->len) - q->pos;
    if (q->pos >= _bnt_annotations->at(repseq).bns->l_pac) {
      continue;
    }
    bns_coor_pac2real(_bnt_annotations->at(repseq).bns, q->pos, j, &seqid);

    // fill in alignment info
    ALIGNMENT refid;
    ParseRefid(_bnt_annotations->at(repseq).bns->anns[seqid].name, &refid);
    refid.strand = q->strand;
    if (q->strand) {
      refid.pos = (refid.start-extend-PAD) + 
	(int)(q->pos - _bnt_annotations->at(repseq).bns->anns[seqid].offset) + 1;
    } else {
      refid.pos = (int)(q->pos - _bnt_annotations->at(repseq).
			bns->anns[seqid].offset - 2*j)+(refid.start-extend-PAD);
    }
    refid.repeat = repseq;
    alignments->push_back(refid);
  }
  return true;
}

bool BWAReadAligner::GetSharedAlns(const vector<ALIGNMENT>& map1,
				  const vector<ALIGNMENT>& map2,
				  vector<ALIGNMENT>* left_refids,
				  vector<ALIGNMENT>* right_refids) {
  if (map1.size() == 0 || map2.size() == 0) return false;
  // get keys that are shared and consistent
  map<int,ALIGNMENT> left_id_to_ref;
  map<int,ALIGNMENT> right_id_to_ref;
  bool found = false;
  for (vector<ALIGNMENT>::const_iterator it = map1.begin();
       it != map1.end(); ++it) {
    // Must see every ID only once
    if (left_id_to_ref.find((*it).id) == left_id_to_ref.end()) {
      left_id_to_ref[(*it).id] = *it;
    } else {
      left_id_to_ref.erase((*it).id);
    }
  }
  for (vector<ALIGNMENT>::const_iterator it = map2.begin();
       it != map2.end(); ++it) {
    // Must see every ID only once
    if (right_id_to_ref.find((*it).id) == right_id_to_ref.end()) {
      right_id_to_ref[(*it).id] = *it;
    } else {
      right_id_to_ref.erase((*it).id);
    }
  }
  for (map<int,ALIGNMENT>::const_iterator it = right_id_to_ref.begin();
       it != right_id_to_ref.end(); ++it) {
    int ref_key = it->first;
    if (left_id_to_ref.find(ref_key) != left_id_to_ref.end()) { // found match!
      // check strand, L, R are compatible
      if (left_id_to_ref[ref_key].strand == it->second.strand &&
	left_id_to_ref[ref_key].left != it->second.left) {
	left_refids->push_back(left_id_to_ref[ref_key]);
	right_refids->push_back(it->second);
	found = true;
      }
    }
  }
  return found;
}

bool BWAReadAligner::AlignMate(const ReadPair& read_pair,
                               vector<ALIGNMENT>* mate_alignments,
                               const string& repseq) {
  const int& num_aligned_read = read_pair.aligned_read_num;
  const string& nucs = reverseComplement(read_pair.reads.at(1-num_aligned_read).orig_nucleotides);
  const string& qual = reverse(read_pair.reads.at(1-num_aligned_read). orig_qual);

  // set up BWA alignment
  int is_comp = _opts->mode&BWA_MODE_COMPREAD;
  bwa_seq_t *seq = (bwa_seq_t*)calloc(1, sizeof(bwa_seq_t));
  seq->bc[0] = 0;
  seq->tid = -1;
  seq->qual = 0;
  seq->full_len = seq->clip_len = seq->len = nucs.length();
  seq->seq = (ubyte_t*)calloc(seq->len, 1);
  seq->qual = (ubyte_t*)calloc(seq->len, 1);
  for (int i = 0; i != seq->full_len; ++i) {
    seq->seq[i] = nst_nt4_table[(int)nucs.at(i)];
    seq->qual[i] = qual.at(i)+33<126?qual.at(i)+33:126;
  }
  seq->rseq = (ubyte_t*)calloc(seq->full_len, 1);
  memcpy(seq->rseq, seq->seq, seq->len);
  seq_reverse(seq->len, seq->seq, 0);
  seq_reverse(seq->len, seq->qual, 0);
  seq_reverse(seq->len, seq->rseq, is_comp);
  seq->name = strdup((const char*)read_pair.
		     reads.at(1-num_aligned_read).ID.c_str());

  // call bwa with appropriate options
  bwa_cal_sa_reg_gap(0, _bwt_references->at(repseq).bwt,
                     1, seq, _default_opts); 

  if (seq->n_aln == 0) {
    bwa_free_read_seq(1, seq);
    if (align_debug) {
      cerr << "[AlignMate]: check mate found 0 alignments " << endl;
    }
    return false;
  }

  // Check alignment coordinates
  if (!GetAlignmentCoordinates(seq, repseq, mate_alignments)) {
    bwa_free_read_seq(1, seq);
    return false;
  }
  bwa_free_read_seq(1, seq);
  return true;
}

bool BWAReadAligner::CheckMateAlignment(const vector<ALIGNMENT>& mate_alignments,
                                        const ALIGNMENT& left_alignment,
                                        const ALIGNMENT& right_alignment,
                                        ALIGNMENT* mate_alignment) {
  // For each, check against STR alignment
  for (vector<ALIGNMENT>::const_iterator it = mate_alignments.begin();
       it != mate_alignments.end(); ++it) {
    // Check chrom, position, and strand
    if (align_debug) {
      cerr << "[CheckMateAlignment]: check mate finds " << it->chrom << " " 
           << (abs(it->pos-left_alignment.pos)) << endl;
    }
    // find start of STR read
    const int& str_pos = left_alignment.left ? 
      left_alignment.pos : right_alignment.pos;
    // check chrom, pos, and strand
    if ((it->chrom == left_alignment.chrom) &&
        (abs(it->pos-str_pos) <= MAX_PAIRED_DIFF) &&
        it->strand != left_alignment.strand) {
      *mate_alignment =  (*it);
      return true;
    }
  }
  return false;
}

bool BWAReadAligner::StitchReads(ReadPair* read_pair,
				 ALIGNMENT* left_alignment,
				 ALIGNMENT* right_alignment) {
  const int& num_aligned_read = read_pair->aligned_read_num;
  const string& seq1 = read_pair->reads.at(num_aligned_read).orig_nucleotides;
  const string& seq2 = reverseComplement(read_pair->reads.at(1-num_aligned_read).
					 orig_nucleotides);
  const string& seq1_qual = read_pair->reads.at(num_aligned_read).orig_qual;
  const string& seq2_qual = reverse(read_pair->reads.
				    at(1-num_aligned_read).orig_qual);
  // Set up 
  vector<float> scores;
  scores.push_back(0);
  float score, max_score = 0;
  size_t overlap_len, max_score_index = -1;

  // Gradually bring ends together and try to stitch
  for (size_t i = 1; i < seq1.length() - MIN_STITCH_OVERLAP; i++) {
    score = 0;
    overlap_len = seq1.length()-i;
    for (size_t j = 0; j < overlap_len; j++) {
      if (j >= seq2.length()) {
	return false;
      }
      if (seq1.at(i+j) == seq2.at(j)) {
	score += 1;
      }
      if (score/overlap_len >= max_score) {
	max_score = score/overlap_len;
	max_score_index = i-1;
      }
    }
    scores.push_back(score/overlap_len);
  }

  // Check if too many matches
  for (size_t i = 0; i < scores.size(); i++) {
    if ((max_score - scores.at(i) <= STITCH_DIFF) && i != max_score_index+1) {
      return false;
    }
  }

  overlap_len = seq1.length() - max_score_index;

  if ((overlap_len < MIN_STITCH_OVERLAP) ||
      (max_score < STITCH_REQUIRED_SCORE)) {
    return false;
  }
  string stitched_string = seq1.substr(0, (int)max_score_index+1);
  string stitched_qual = seq1_qual.substr(0, (int)max_score_index+1);
  string na, nb, qa, qb;
  for (size_t i = 0; i < overlap_len; i++) {
    na = seq1.substr(max_score_index+i+1,1);
    nb = seq2.substr(i,1);
    qa = seq1_qual.substr(max_score_index+i+1,1);
    qb = seq2_qual.substr(i,1);
    if (qa > qb) {
      stitched_string.append(na);
      stitched_qual.append(qa);
    } else if (qa < qb) {
      stitched_string.append(nb);
      stitched_qual.append(qb);
    } else {
      stitched_string.append(na);
      stitched_qual.append(qa);
    }
  }
  stitched_string.append(seq2.substr(overlap_len));
  stitched_qual.append(seq2_qual.substr(overlap_len));
  
  // put stitched info in aligned read
  read_pair->reads.at(num_aligned_read).nucleotides = stitched_string;
  read_pair->reads.at(num_aligned_read).quality_scores = stitched_qual;
  read_pair->reads.at(num_aligned_read).right_flank_nuc = 
    stitched_string.substr(seq1.length()-read_pair->reads.at(num_aligned_read).
			   right_flank_index_from_end -
			   read_pair->reads.at(num_aligned_read).
			   right_flank_nuc.length());
  if (!left_alignment->left) {
    right_alignment->pos -= ((seq2.length() - overlap_len)+
			     read_pair->reads.at(num_aligned_read).
			     right_flank_index_from_end);
  } else {
    left_alignment->pos -= read_pair->reads.at(num_aligned_read).
      left_flank_index_from_start;
  }
  return true;
}

bool BWAReadAligner::OutputAlignment(ReadPair* read_pair,
				     const ALIGNMENT& left_alignment,
				     const ALIGNMENT& right_alignment,
				     const ALIGNMENT& mate_alignment,
				     bool treat_as_paired) {
  const int& aligned_read_num = read_pair->aligned_read_num;
  read_pair->treat_as_paired = treat_as_paired;
  // Set info for aligned read
  read_pair->reads.at(aligned_read_num).chrom = left_alignment.chrom;
  read_pair->reads.at(aligned_read_num).strid = left_alignment.id;
  if (left_alignment.left) {
    read_pair->reads.at(aligned_read_num).msStart = left_alignment.start;
    read_pair->reads.at(aligned_read_num).msEnd = left_alignment.end;
  } else {
    read_pair->reads.at(aligned_read_num).msStart = right_alignment.start;
    read_pair->reads.at(aligned_read_num).msEnd = right_alignment.end;
  }
  read_pair->reads.at(aligned_read_num).refCopyNum = left_alignment.copynum;
  read_pair->reads.at(aligned_read_num).reverse = !left_alignment.left;

  read_pair->reads.at(aligned_read_num).lStart = left_alignment.pos;
  read_pair->reads.at(aligned_read_num).lEnd = left_alignment.pos +
    read_pair->reads.at(aligned_read_num).left_flank_nuc.length();
  read_pair->reads.at(aligned_read_num).rStart = right_alignment.pos;
  read_pair->reads.at(aligned_read_num).rEnd = right_alignment.pos +
    read_pair->reads.at(aligned_read_num).right_flank_nuc.length();

  // checks to make sure the alignment is reasonable
  // coords make sense on forward strand
  if (!read_pair->reads.at(aligned_read_num).reverse & 
      !read_pair->reads.at(aligned_read_num).partial 
      & ((read_pair->reads.at(aligned_read_num).lStart >=
	  read_pair->reads.at(aligned_read_num).msStart) 
	 | (read_pair->reads.at(aligned_read_num).rEnd <=
	    read_pair->reads.at(aligned_read_num).msEnd))) {
    return false;
  }
  // coords make sense on reverse strand
  if (read_pair->reads.at(aligned_read_num).reverse & 
      !read_pair->reads.at(aligned_read_num).partial 
      & ((read_pair->reads.at(aligned_read_num).rStart >=
	  read_pair->reads.at(aligned_read_num).msStart)
	 | (read_pair->reads.at(aligned_read_num).lEnd <=
	    read_pair->reads.at(aligned_read_num).msEnd))) {
    return false;
  }

  // coords of flanks are positive
  if ((read_pair->reads.at(aligned_read_num).rStart <0) |
      (read_pair->reads.at(aligned_read_num).lStart <0) | 
      (read_pair->reads.at(aligned_read_num).rEnd <0) |
      (read_pair->reads.at(aligned_read_num).lEnd < 0) | 
      (read_pair->reads.at(aligned_read_num).msStart < 0)) {
    return false;
  }

  // Adjust alignment and STR call
  try {
    if (!AdjustAlignment(&read_pair->reads.at(aligned_read_num),
			 read_pair->reads.at(aligned_read_num).partial,
			 !read_pair->reads.at(aligned_read_num).left_all_repeats,
			 !read_pair->reads.at(aligned_read_num).right_all_repeats)) {
      return false;
    }
  } catch (std::out_of_range & exception) {
    return false;
  }

  // Make sure alignment meets requirements
  if (unit && !read_pair->reads.at(aligned_read_num).partial) {
    if (read_pair->reads.at(aligned_read_num).diffFromRef %
	read_pair->reads.at(aligned_read_num).ms_repeat_best_period != 0)
      return false;
  }
  if (((abs(read_pair->reads.at(aligned_read_num).diffFromRef) > max_diff_ref) &&
	 !read_pair->reads.at(aligned_read_num).partial)) {
    return false;
  }

  if (treat_as_paired) {
    read_pair->reads.at(1-aligned_read_num).read_start = mate_alignment.pos;
    read_pair->reads.at(1-aligned_read_num).sw_score = 255; // TODO
    read_pair->reads.at(1-aligned_read_num).reverse =
      !read_pair->reads.at(aligned_read_num).reverse;
    // get cigar
    CIGAR_LIST cigar_list;
    string aln_seq, ref_seq;
    int score;
    const size_t& reglen = read_pair->reads.at(1-aligned_read_num).nucleotides.length();
    const REFSEQ& refseq = _ref_sequences->at(read_pair->reads.at(aligned_read_num).strid);
    const size_t& start_pos = read_pair->reads.at(1-aligned_read_num).reverse ?
      mate_alignment.pos : mate_alignment.pos-1;
    const string& rseq = refseq.sequence.substr(start_pos - refseq.start, reglen);
    const string& aseq = read_pair->reads.at(1-aligned_read_num).reverse ?
      reverseComplement(read_pair->reads.at(1-aligned_read_num).nucleotides) :
      read_pair->reads.at(1-aligned_read_num).nucleotides;
    nw(aseq, rseq, aln_seq, ref_seq, false, &score, &cigar_list);

    read_pair->reads.at(1-aligned_read_num).cigar_string =
      cigar_list.cigar_string;
    read_pair->reads.at(1-aligned_read_num).cigar =
      cigar_list.cigars;
  }
  return true;
}

bool BWAReadAligner::AdjustAlignment(MSReadRecord* aligned_read, bool partial,
				     bool left_aligned, bool right_aligned) {
  // get reference sequence
  const size_t& reglen = !aligned_read->reverse ? 
    (aligned_read->rEnd - aligned_read->lStart) :
    (aligned_read->lEnd - aligned_read->rStart);
  if (_ref_sequences->find(aligned_read->strid) == 
      _ref_sequences->end()) return false;
  const REFSEQ& refseq = _ref_sequences->at(aligned_read->strid);
  size_t start_pos = !aligned_read->reverse ? 
    aligned_read->lStart-1 : aligned_read->rStart;
  // start is in rep region, use end other end
  if (partial) {
    if (aligned_read->reverse & left_aligned) {
      start_pos = aligned_read->lEnd - reglen-1;
    } 
    if (!aligned_read->reverse & right_aligned) {
      start_pos = aligned_read->rEnd - reglen-1;
    }
  }

  // get sequences to pairwise align
  const string& rseq = refseq.sequence.substr(start_pos - refseq.start, reglen);
  const string& aligned_seq = !aligned_read->reverse ? 
    aligned_read->nucleotides : reverseComplement(aligned_read->nucleotides);

  // update coords
  aligned_read->read_start = start_pos;
  aligned_read->read_end = start_pos + reglen;

  if (debug_adjust) {
    cerr << "[AdjustAlignment]: Partial " << partial 
         << " Reverse " << aligned_read->reverse 
         << " left rep " << !left_aligned 
         << " right rep " << !right_aligned << endl;
    cerr << "[AdjustAlignment]: " << rseq << endl 
         << "[AdjustAlignment]: " << aligned_seq << endl;
  }

  // Global alignment read vs. region aligned to 
  string aligned_seq_sw, ref_seq_sw;
  int sw_score;
  CIGAR_LIST cigar_list;
  nw(aligned_seq, rseq, aligned_seq_sw, ref_seq_sw,
     false, &sw_score, &cigar_list);

  // Fix off by one alignment problem
  if (cigar_list.cigars.at(0).num == 1) {
    if (cigar_list.cigars.at(0).cigar_type=='I') {
      cigar_list.cigars.erase(cigar_list.cigars.begin());
      cigar_list.cigars.at(1).num++;
      aligned_read->read_start++;
    }
  }

  // get rid of end gaps
  if (cigar_list.cigars.at(0).cigar_type == 'D') {
    const int& num = cigar_list.cigars.at(0).num;
    aligned_read->read_start += num;
    cigar_list.cigars.erase(cigar_list.cigars.begin());
    cigar_list.ResetString();
  }
  if (cigar_list.cigars.at(cigar_list.cigars.size() -1).cigar_type == 'D') {
    cigar_list.cigars.erase(cigar_list.cigars.end()-1);
    cigar_list.ResetString();
  }
  
  // Update info in aligned read
  aligned_read->sw_score = sw_score; // TODO change
  aligned_read->cigar = cigar_list.cigars;
  aligned_read->cigar_string = cigar_list.cigar_string;

  if (debug_adjust) {
    cerr << "[AdjustAlignment]: " << aligned_read->ID << endl;
    cerr << "[AdjustAlignment]: " << aligned_seq_sw << endl;
    cerr << "[AdjustAlignment]: " << ref_seq_sw << endl;
    cerr << "[AdjustAlignment]: " << sw_score << endl;
    cerr << "[AdjustAlignment]: " << cigar_list.cigar_string << endl;
    cerr << "[AdjustAlignment]: " << aligned_read->diffFromRef << endl;
  }

  // Check if alignment is reasonably good
  if (sw_score < min_sw_score ) return false;

  // Readjust if partial
  if (partial) {
    if (!AdjustPartialAlignment(aligned_read,cigar_list,
				left_aligned, right_aligned,
				start_pos, reglen)) {
      return false;
    }
  } else {
    // Update STR allele
    if (!GetSTRAllele(aligned_read, cigar_list)) return false;
  }
  return true;
}

bool BWAReadAligner::AdjustPartialAlignment(MSReadRecord* aligned_read,
					    const CIGAR_LIST& cigar_list,
					    bool left_aligned, bool right_aligned,
					    int start_pos, int reglen) {
  // how many base pairs are spanned in the reference by this read?
  int span = 0;
  // Index into cigar string
  size_t i;
  // how many base pairs are part of the actual str sequence?
  int strbp;
  // Difference in bp from the reference
  int diff_from_ref = 0;
  // Length of the STR region
  int ms_length = aligned_read->msEnd - aligned_read->msStart;

  // Get span of STR region
  for (i = 0; i < cigar_list.cigars.size(); i++) {
    if ((cigar_list.cigars.at(i).cigar_type == 'M') |
	(cigar_list.cigars.at(i).cigar_type == 'D')) {
      span += cigar_list.cigars.at(i).num;
    }
    if (cigar_list.cigars.at(i).cigar_type == 'D') {
      diff_from_ref -= cigar_list.cigars.at(i).num;
    }
    if (cigar_list.cigars.at(i).cigar_type == 'I') {
      diff_from_ref += cigar_list.cigars.at(i).num;
    }
  }

  // Check if actually partially spanning
  // Case 1: starts at left. see if total span is > (str index+ms len)
  if ((!aligned_read->reverse & left_aligned) |
      (aligned_read->reverse & right_aligned)) {
    strbp = span - (aligned_read->msStart - aligned_read->read_start);
    // Check if end extends beyond the read
    if (aligned_read->msEnd - aligned_read->read_start >= span) {
      // Check if reasonable
      if (strbp+diff_from_ref < 0 || 
	  strbp+diff_from_ref >= (int)aligned_read->nucleotides.length() ||
	  strbp < (int)MIN_STR_LENGTH) return false;
      // Update read data
      aligned_read->partial = true;
      aligned_read->diffFromRef = (strbp+diff_from_ref) - ms_length ;
      aligned_read->detected_ms_nuc = aligned_read->reverse ? 
	aligned_read->nucleotides.substr(0,strbp+diff_from_ref) : 
	reverseComplement(aligned_read->nucleotides).
	substr(0,strbp+diff_from_ref);	
      return true;
    } else {
      if (debug_adjust) {
        cerr << "[AdjustPartialAlignment]: actually not partial" << endl;
      }
      aligned_read->partial = false;
      return true;
    }
  }

  // Case 2: anchored at the right
  if ((!aligned_read->reverse & right_aligned) |
      (aligned_read->reverse & left_aligned)) {
    strbp = span - (start_pos+reglen-aligned_read->msEnd);
    // Check if beginning extends beyond the msStart
    if ((start_pos  + reglen - aligned_read->msStart) >= span) {
      // check if reasonable
      if (strbp+diff_from_ref < 0 || 
	  strbp+diff_from_ref >= (int)aligned_read->nucleotides.length() ||
	  strbp <= (int)MIN_STR_LENGTH) {
	return false;
      }
      // update read info
      aligned_read->partial = true;
      aligned_read->diffFromRef = (strbp +diff_from_ref)-ms_length;
      aligned_read->detected_ms_nuc = aligned_read->reverse ? 
	aligned_read->nucleotides.
	substr(aligned_read->nucleotides.length() -
	       strbp-diff_from_ref, strbp+diff_from_ref) :
	reverseComplement(aligned_read->nucleotides).
	substr(aligned_read->nucleotides.length() -
	       strbp-diff_from_ref,strbp+diff_from_ref);	
      return true;
    } else {
      if (debug_adjust) {
        cerr << "[AdjustPartialAlignment]: actually not partial" << endl;
      }
      aligned_read->partial = false;
      return true;
    }
  }
  return false;
}

bool BWAReadAligner::GetSTRAllele(MSReadRecord* aligned_read,
				  const CIGAR_LIST& cigar_list) {
  const bool& partial = aligned_read->partial;
  // index where STR starts in the read
  size_t str_index = aligned_read->msStart-aligned_read->read_start +1;
  // Length of the total STR region
  size_t ms_length = aligned_read->msEnd - aligned_read->msStart -1;
  // If right end anchored
  if (partial & (str_index < 0)) {
    str_index = 0;
    ms_length = aligned_read->msEnd - aligned_read->read_start;
  }
  // If left end anchored
  if (partial & (str_index+ms_length >= aligned_read->nucleotides.length())) {
    ms_length = aligned_read->nucleotides.length() - str_index - 1;
  }

  // If alignment is too messy, get rid of it
  if (cigar_list.cigars.size() > MAX_CIGAR_SIZE) {
    return false;
  }

  // same as reference
  if (cigar_list.cigars.size() == 1) {
    if (aligned_read->reverse) {
      aligned_read->detected_ms_nuc =
	reverseComplement(aligned_read->nucleotides).
	substr(str_index, ms_length);
    } else {
      aligned_read->detected_ms_nuc =
	aligned_read->nucleotides.substr(str_index, ms_length);
    }
    aligned_read->diffFromRef = 0;
    return (aligned_read->detected_ms_nuc.length() > MIN_STR_LENGTH);
  }

  // get only cigar score spanning the STR
  const int& str_start_in_cigar = aligned_read->msStart - aligned_read->read_start + 1;
  int pos = 0;
  int bp = 0;
  size_t cigar_index = 0;
  int diff = 0;
  int total_cigar_pos = 0;
  int diff_from_ref = 0;
  char cigar_type;
  CIGAR_LIST new_cigar_list;
  CIGAR_LIST str_cigar_list;
 
  // get rid of left flanking region
  while (pos <= str_start_in_cigar && 
	 cigar_index < cigar_list.cigars.size()) {
    bp = cigar_list.cigars.at(cigar_index).num;
    cigar_type = cigar_list.cigars.at(cigar_index).cigar_type;
    // If match or del, increment position
    if (cigar_type == 'M' || cigar_type == 'D') pos += bp;
    // bp to go until we hit STR
    diff = pos - str_start_in_cigar + 1;
    if (diff >= 0) {
      if (cigar_index == 0) {
	new_cigar_list.cigars.resize(cigar_list.cigars.size() - cigar_index + 1);
	copy(cigar_list.cigars.begin()+cigar_index, cigar_list.cigars.end(),
	     new_cigar_list.cigars.begin());
	diff -= cigar_list.cigars.at(cigar_index).num;
      } else {
	new_cigar_list.cigars.resize(cigar_list.cigars.size() - cigar_index + 2);
	copy(cigar_list.cigars.begin()+cigar_index-1, cigar_list.cigars.end(),
	     new_cigar_list.cigars.begin());
	diff -= cigar_list.cigars.at(cigar_index-1).num;
      }
      break;
    }
    cigar_index += 1;
  }

  // Update STR cigar taking away left flank
  str_cigar_list.cigars = new_cigar_list.cigars;
  str_cigar_list.ResetString();
  new_cigar_list.cigars.clear();

  // get rid of right flank cigars
  total_cigar_pos = ms_length;
  cigar_index = 0;
  pos = diff;
  while(pos < total_cigar_pos) {
    if (cigar_index >= str_cigar_list.cigars.size()) {
      return false;
    }
    bp = str_cigar_list.cigars.at(cigar_index).num;
    cigar_type = str_cigar_list.cigars.at(cigar_index).cigar_type;
    if (cigar_type == 'M' || cigar_type == 'D') pos += bp;
    diff = pos-total_cigar_pos;
    if (diff >= 0) {
      new_cigar_list.cigars.resize(cigar_index+1);
      copy(str_cigar_list.cigars.begin(), str_cigar_list.cigars.begin()+cigar_index+1,
	   new_cigar_list.cigars.begin());
      break;
    }
    cigar_index += 1;
  }
  str_cigar_list.cigars.clear();
  str_cigar_list.cigars = new_cigar_list.cigars;

  // set diff from ref
  diff_from_ref = 0;
  for (size_t i = 0; i < str_cigar_list.cigars.size(); i++) {
    if (str_cigar_list.cigars.at(i).cigar_type == 'I') {
      diff_from_ref += str_cigar_list.cigars.at(i).num;
    }
    if (str_cigar_list.cigars.at(i).cigar_type == 'D') {
      diff_from_ref -= str_cigar_list.cigars.at(i).num;
    }
  }

  // set STR region
  string ms_nuc;
  if (aligned_read->reverse) {
    string rev_read = reverseComplement(aligned_read->nucleotides);
    ms_nuc =  rev_read.substr(str_index, ms_length+diff_from_ref);
  } else {
    ms_nuc =  aligned_read->nucleotides.substr(str_index, ms_length+diff_from_ref);
  }
  if (ms_nuc.length() <= MIN_STR_LENGTH) return false;
  aligned_read->diffFromRef = diff_from_ref;
  aligned_read->detected_ms_nuc = ms_nuc;
  return true;
}


BWAReadAligner::~BWAReadAligner(){}

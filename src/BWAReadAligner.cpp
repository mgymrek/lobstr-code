/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdexcept>

#include "BWAReadAligner.h"
#include "cigar.h"
#include "common.h"
#include "nw.h"
#include "runtime_parameters.h"


extern unsigned char nst_nt4_table[256];
int pad=50;
int extend_partial = 0;
size_t MAX_CIGAR_SIZE = 6;
size_t MAX_ALLOWED_FLANK = 10;

static size_t count(string s, char c) {
  size_t num = 0;
  for (size_t i = 0; i < s.length(); i++) {
    if (s.at(i) == c) num++;
  }
  return num;
}

BWAReadAligner::BWAReadAligner(map<std::string, BWT>* bwt_references,
			       map<std::string, BNT>* bnt_annotations,
			       map<int, REFSEQ>* ref_sequences,
			       gap_opt_t *opts) {
  bwase_initialize();
  _bwt_references = bwt_references;
  _bnt_annotations = bnt_annotations;
  _ref_sequences = ref_sequences;
  _opts = opts;
}

bool BWAReadAligner::ProcessRead(MSReadRecord* read) {
  bool left_all_repeats = false;
  bool right_all_repeats = false;
  bool partial = false;
  read->partial = false;
  if (align_debug)
    cout << endl << "processing " << read->ID << " " << read->nucleotides << " repseq " << read->repseq << endl;
  const string& repseq = read->repseq;

  // check that we have that repseq
  if (_bwt_references->find(repseq) == _bwt_references->end()) {
    return false;
  }

  // check that the repseq is valid
  if (count(repseq,repseq.at(0)) == repseq.length()) return false;
  
  // check that it's not a perect repeat
  if (IsPerfectRepeat(read->nucleotides, read->repseq)) {
    return false;
  }

  // Check if flanks are perfect repeats
  if (IsPerfectRepeat(read->left_flank_nuc, read->repseq)) 
    left_all_repeats = true;
  if (IsPerfectRepeat(read->right_flank_nuc, read->repseq))
    right_all_repeats = true;
  // don't continue if both are fully repetitive
  if (left_all_repeats && right_all_repeats) return false;

  // set up bwa
  int is_comp = _opts->mode&BWA_MODE_COMPREAD;

  // get bwa_seq_t *seqs
  bwa_seq_t *seqs, *seq_left, *seq_right;
  seqs=(bwa_seq_t*)calloc(2, sizeof(bwa_seq_t));
  seq_left = &seqs[0];
  seq_right= &seqs[1];

  // left in for historical reasons...
  read->left_flank_nuc = reverseComplement(read->left_flank_nuc);
  read->right_flank_nuc = reverseComplement(read->right_flank_nuc);

  // LEFT FLANKING REGION
  const char* left_flank = read->left_flank_nuc.c_str();
  const char* left_qual =  read->quality_scores.substr(0, strlen(left_flank)).c_str();

  seq_left->bc[0] = 0;
  seq_left->tid = -1;
  seq_left->qual = 0;
  seq_left->full_len = seq_left->clip_len = seq_left->len = strlen(left_flank);
  seq_left->seq = (ubyte_t*)calloc(seq_left->len, 1);
  seq_left->qual = (ubyte_t*)calloc(seq_left->len, 1);
  for (int i = 0; i != seq_left->full_len; ++i) {
    seq_left->seq[i] = nst_nt4_table[(int)left_flank[i]];
    if (fastq) {
      seq_left->qual[i] = left_qual[i]+33<126?left_qual[i]+33:126;
    } else {
      seq_left->qual[i] = 126;
    }
  }
  seq_left->rseq = (ubyte_t*)calloc(seq_left->full_len,1);
  memcpy(seq_left->rseq, seq_left->seq, seq_left->len);
  seq_reverse(seq_left->len, seq_left->seq, 0);
  seq_reverse(seq_left->len, seq_left->rseq, is_comp);
  seq_left->name = strdup((const char*)read->ID.c_str());
				   
  // RIGHT FLANKING REGION
  const char* right_flank = read->right_flank_nuc.c_str();
  const char* right_qual = read->quality_scores.substr(read->quality_scores.size() -
						       strlen(right_flank),
						       strlen(right_flank)).c_str();
  seq_right->bc[0] = 0;
  seq_right->tid = -1;
  seq_right->qual = 0;
  seq_right->full_len = seq_right->clip_len = seq_right->len = strlen(right_flank);
  seq_right->seq = (ubyte_t*)calloc(seq_right->len, 1);
  seq_right->qual = (ubyte_t*)calloc(seq_right->len, 1);
  for (int i = 0; i != seq_right->full_len; ++i) {
    seq_right->seq[i] = nst_nt4_table[(int)right_flank[i]];
    if (fastq) {
      seq_right->qual[i] = right_qual[i]+33<126?right_qual[i]+33:126;
    } else {
      seq_right->qual[i] = 126;
    }
  }
  seq_right->rseq = (ubyte_t*)calloc(seq_right->full_len,1);
  memcpy(seq_right->rseq, seq_right->seq, seq_right->len);
  seq_reverse(seq_right->len, seq_right->seq, 0);
  seq_reverse(seq_right->len, seq_right->rseq, is_comp);
  seq_right->name = strdup((const char*)read->ID.c_str());

  if (align_debug)
    cout << "left flank " << left_flank << " right flank " << right_flank << endl;

  // call bwa with appropriate options
  bwa_cal_sa_reg_gap(0, _bwt_references->at(repseq).bwt, 2, seqs, _opts); 

  // Find matches between them (need same STR, same direction)
  int num_left_aln = seq_left->n_aln;
  int num_right_aln = seq_right->n_aln;
  if (num_left_aln == 0 || num_right_aln == 0) {
    bwa_free_read_seq(2, seqs);
    return false;
  }

  // get ref seqs and index of alignment for each
  list<ALIGNMENT> left_alignments;
  list<ALIGNMENT> right_alignments;

  // look for left alignment
  if (!left_all_repeats) {
    if (!GetAlignmentCoordinates(seq_left, repseq, &left_alignments)) {
      bwa_free_read_seq(2, seqs);
      return false;
    }
  }

  // look for right alignment
  if (!right_all_repeats) {
    if (!GetAlignmentCoordinates(seq_right, repseq, &right_alignments)) {
      bwa_free_read_seq(2, seqs);
      return false;
    }
  }

  // get shared keys
  ALIGNMENT left_refid;
  ALIGNMENT right_refid;
  bool left_cleared = false;
  bool right_cleared = false;
  if (!GetSharedAln(left_alignments, right_alignments,
  		    &left_refid, &right_refid)) {
    // if one flank maps everywhere, clear it
    if (left_alignments.size() > MAX_ALLOWED_FLANK) {
      left_alignments.clear();
    }
    if (right_alignments.size() > MAX_ALLOWED_FLANK) {
      right_alignments.clear();
    }
    // check if one flanking region aligns
    if (left_alignments.size() + right_alignments.size() == 1) {
      partial = true;
      read->partial = true;
      // Case 1: anchored the left side
      if (left_alignments.size() != 0 && (right_all_repeats||right_cleared)) {
	right_all_repeats = true;
	left_refid = left_alignments.front();
	right_refid.left = !left_refid.left;
	right_refid.start = left_refid.start - extend;
	right_refid.end = left_refid.end - extend;
	// set right pos
	if (!left_refid.left) {
	  // reverse complement (rev 1 right rep 1)
	  right_refid.pos = left_refid.pos + strlen(left_flank) - read->nucleotides.size();
	} else {
	  // rev 0 right rep 1
	  right_refid.pos = left_refid.pos + read->nucleotides.size() - strlen(right_flank);
	}
	if (debug_adjust) {
	  cout << "Partial, right rep " << read->nucleotides << " " << " reverse " << !left_refid.left << " left rep " << 0 << " right rep " << 1 << endl;
	}
      }
      // Case 2: anchored the right side
      else if (right_alignments.size() != 0 && (left_all_repeats||left_cleared)) {
	left_all_repeats = true;
	right_refid = right_alignments.front();
	left_refid.left = !right_refid.left;
	left_refid.chrom = right_refid.chrom;
	left_refid.id = right_refid.id;
	left_refid.start = right_refid.start - extend;
	left_refid.end = right_refid.end - extend;
	left_refid.copynum = right_refid.copynum;
	// set left pos
	if (!left_refid.left) {
	  // reverse complement (rev 1 left 1)
	  left_refid.pos = right_refid.pos + read->nucleotides.size() - strlen(left_flank);
	} else { // rev 0 right 1
	  left_refid.pos = right_refid.pos + strlen(right_flank) - read->nucleotides.size();
	}
	if (debug_adjust) {
	  cout << "Partial, left rep " << read->nucleotides << " " << " reverse " << !left_refid.left << " left rep " << 1 << " right rep " << 0 << endl;
	}
      } else {
	bwa_free_read_seq(2, seqs);
	return false;
      }
    } else {
      bwa_free_read_seq(2, seqs);
      return false;
    }
  }
  if (align_debug)
    cout << "found shared aln " <<endl;


  // get the desired info out of seqs and put it back into the MSReadRecord
  read->chrom = left_refid.chrom;
  read->strid = left_refid.id;
  if (left_refid.left) {
    read->msStart = left_refid.start;
    read->msEnd = left_refid.end;
  } else {
    read->msStart = right_refid.start;
    read->msEnd = right_refid.end;
  }
  read->msRepeat = repseq;
  read->refCopyNum = left_refid.copynum;
  read->reverse = !left_refid.left;
  read->lStart = left_refid.pos;
  read->lEnd = left_refid.pos + strlen(left_flank);
  read->rStart = right_refid.pos;
  read->rEnd = right_refid.pos + strlen(right_flank);
  bwa_free_read_seq(2, seqs);

  // complement back again
  read->left_flank_nuc = reverseComplement(read->left_flank_nuc);
  read->right_flank_nuc = reverseComplement(read->right_flank_nuc);


  // checks to make sure the alignment is reasonable
  if (!read->reverse & !partial & ((read->lStart >= read->msStart) | (read->rEnd <= read->msEnd))) {
    return false;
  }
  if (read->reverse & !partial & ((read->rStart >= read->msStart) | (read->lEnd <= read->msEnd))) {
    return false;
  }
  if ((read->rStart <0) | (read->lStart <0) | (read->rEnd <0) | (read->lEnd < 0) | (read->msStart < 0)) {
    return false;
  }

  // adjust alignment and set cigar score
  if (adjust) {
    try {
      if (!AdjustAlignment(read, partial, !left_all_repeats, !right_all_repeats)) {
	return false;
      }
    } catch (std::out_of_range & exception) {
      return false;
    }
  }
  if (unit) {
    if (read->diffFromRef % read->ms_repeat_best_period != 0)
      return false;
  }
  return ((abs(read->diffFromRef) <= max_diff_ref)||partial);
}

bool BWAReadAligner::GetSharedAln(const list<ALIGNMENT>& map1,
				  const list<ALIGNMENT>& map2,
				  ALIGNMENT* left_refid,
				  ALIGNMENT* right_refid) {
  if (map1.size() == 0 || map2.size() == 0) return false;
  // get keys that are shared and consistent
  map<int,ALIGNMENT> left_id_to_ref;
  bool found = false;
  for (list<ALIGNMENT>::const_iterator it = map1.begin();
       it != map1.end(); ++it) {
    if (left_id_to_ref.find((*it).id) == left_id_to_ref.end()) {
      left_id_to_ref[(*it).id] = *it;
    } else {
      left_id_to_ref.erase((*it).id);
    }
  }
  for (list<ALIGNMENT>::const_iterator it = map2.begin();
       it != map2.end(); ++it) {
    int ref_key = (*it).id;
    if (left_id_to_ref.find(ref_key) != left_id_to_ref.end()) { // found match!
      // check strand, L, R are compatible
      /* Need:
	 1. Strand the same
	 2. One right and one left
	 3. If L, R -> repeat is in forward dir. 
	    If R, L -> repeat is in backward dir 
	    (check #3 in calling function) TODO
       */
      if (left_id_to_ref[ref_key].strand == (*it).strand &&
	left_id_to_ref[ref_key].left != (*it).left) {
	if (found) return false; // multiple
	found = true;
	*left_refid = left_id_to_ref[ref_key];
	*right_refid = *it;
      }
    }
  }
  return found;
 }

void BWAReadAligner::ParseRefid(const string& refstring, ALIGNMENT* refid) {
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
}

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
					     list<ALIGNMENT>* alignments) {
  // fill in alignment properties
  bwa_aln2seq_core(aligned_seqs->n_aln, aligned_seqs->aln, aligned_seqs,0,1000);

  if (align_debug)
    cout << "nmulti " << aligned_seqs->n_multi << endl;

  // forward
  if (aligned_seqs->strand) {
    bwa_cal_pac_pos_core(_bwt_references->at(repseq).bwt[0], 0, &aligned_seqs[0],
			 _opts->max_diff, _opts->fnr);
  }
  for (int j = 0; j < aligned_seqs[0].n_multi; ++j) {
    bwt_multi1_t *p = aligned_seqs[0].multi + j;
    if (p->strand) p->pos = bwt_sa(_bwt_references->at(repseq).bwt[0], p->pos);
  }

  // reverse
  if (!aligned_seqs[0].strand) {
    bwa_cal_pac_pos_core(0, _bwt_references->at(repseq).bwt[1], &aligned_seqs[0], 
			 _opts->max_diff, _opts->fnr);
  }
  for (int j = 0; j < aligned_seqs[0].n_multi; ++j) {
    bwt_multi1_t *p = aligned_seqs[0].multi + j;
    if (!p->strand) p->pos = _bwt_references->at(repseq).bwt[1]->seq_len -
		      bwt_sa(_bwt_references->at(repseq).bwt[1], p->pos) +
		      aligned_seqs[0].len;
  }

  // get coords
  for (int i = 0; i < aligned_seqs->n_multi; ++i) {
    bwt_multi1_t *q = aligned_seqs->multi + i;
    int j, nn, seqid; 
    j = pos_end_multi(q, aligned_seqs->len) - q->pos;
    if (q->pos >= _bnt_annotations->at(repseq).bns->l_pac) {
      continue;
    }
    nn = bns_coor_pac2real(_bnt_annotations->at(repseq).bns, q->pos, j, &seqid);
    ALIGNMENT refid;
    ParseRefid(_bnt_annotations->at(repseq).bns->anns[seqid].name, &refid);

    if (align_debug) {
      cout << "adding " << _bnt_annotations->at(repseq).bns->anns[seqid].name << endl;
    }    
    refid.strand = q->strand;
    if (align_debug) {
      cout << "start " << refid.start << " qpos " << q->pos << " annotation " 
	   << _bnt_annotations->at(repseq).bns->anns[seqid].offset << endl;
    }
    if (q->strand) {
      refid.pos = (refid.start-extend-pad) + 
	(int)(q->pos - _bnt_annotations->at(repseq).bns->anns[seqid].offset) + 1;
    } else {
      refid.pos = (int)(q->pos - _bnt_annotations->at(repseq).bns->anns[seqid].offset - 2*j)
	+ (refid.start-extend-pad);
    }
    alignments->push_back(refid);
  }
  return true;
}

/*
Perform local realignment
Adjust exact STR boundaries
Update cigar score
 */
bool BWAReadAligner::AdjustAlignment(MSReadRecord* aligned_read, bool partial, bool left_aligned, bool right_aligned) {
  // get reference sequence
  size_t reglen = !aligned_read->reverse ? (aligned_read->rEnd - aligned_read->lStart) :
    (aligned_read->lEnd - aligned_read->rStart);
  if (_ref_sequences->find(aligned_read->strid) == 
      _ref_sequences->end()) return false;
  REFSEQ refseq = _ref_sequences->at(aligned_read->strid);
  size_t start_pos = !aligned_read->reverse ? aligned_read->lStart-1 : aligned_read->rStart;
  if (start_pos-refseq.start < 0 || 
      reglen + (start_pos - refseq.start) >= refseq.sequence.length()) {return false;}
  // start is in rep region, use end other end
  if (partial) {
    if (aligned_read->reverse & left_aligned) {
      start_pos = aligned_read->lEnd - reglen-1;
    } 
    if (!aligned_read->reverse & right_aligned) {
      start_pos = aligned_read->rEnd - reglen-1;
    }
  }
  const string& rseq = refseq.sequence.substr(start_pos - refseq.start, reglen);
  const string& aligned_seq = !aligned_read->reverse ? aligned_read->nucleotides : reverseComplement(aligned_read->nucleotides);

  // update coords
  aligned_read->read_start = start_pos;
  aligned_read->read_end = start_pos + reglen;
  if (aligned_read->reverse) {
    aligned_read->orig_start = start_pos - aligned_read->right_flank_index_from_end;
  } else {
    aligned_read->orig_start = start_pos - aligned_read->left_flank_index_from_start;
  }
  if (debug_adjust) {
    cout << endl << "Partial " << partial << " Reverse " << aligned_read->reverse << " left rep " << !left_aligned << " right rep " << !right_aligned << endl;
    cout << rseq << " " << aligned_seq << endl;
  }

  // Global alignment read vs. region aligned to 
  string aligned_seq_sw;
  string ref_seq_sw;
  int sw_score;
  CIGAR_LIST cigar_list;
  nw(aligned_seq, rseq, aligned_seq_sw, ref_seq_sw, false, &sw_score, &cigar_list);
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
    int num = cigar_list.cigars.at(0).num;
    aligned_read->read_start += num;
    cigar_list.cigars.erase(cigar_list.cigars.begin());
    cigar_list.ResetString();
  }
  if (cigar_list.cigars.at(cigar_list.cigars.size() -1).cigar_type == 'D') {
    cigar_list.cigars.erase(cigar_list.cigars.end()-1);
    cigar_list.ResetString();
  }
  

  aligned_read->sw_score = sw_score;
  aligned_read->cigar = cigar_list.cigars;
  aligned_read->cigar_string = cigar_list.cigar_string;

  if (debug_adjust) {
    cout << aligned_read->ID << endl;
    cout << aligned_seq_sw << endl;
    cout << ref_seq_sw << endl;
    cout << sw_score << endl;
    cout << cigar_list.cigar_string << endl;
    cout << aligned_read->diffFromRef << endl;
  }

  // Check if alignment is reasonably good
  if (sw_score < min_sw_score ) return false;

  // Adjust if partial
  if (partial) {
    if (!AdjustPartialAlignment(aligned_read,cigar_list, left_aligned, right_aligned, start_pos, reglen)) {
      return false;
    }
  }

  if (!GetSTRAllele(aligned_read, cigar_list, &partial)) return false;
  if (debug_adjust) {
    cout << aligned_read->diffFromRef << " " << aligned_read->detected_ms_nuc << " "
	 << aligned_read->reverse << endl << endl;
  }

  return true;
}

bool BWAReadAligner::AdjustPartialAlignment(MSReadRecord* aligned_read, const CIGAR_LIST& cigar_list, bool left_aligned, bool right_aligned, int start_pos, int reglen) {

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
  // Chack if actually partially spanning
  // Case 1: starts at left. see if total span is > (str index+ms len)
  if ((!aligned_read->reverse & left_aligned) | (aligned_read->reverse & right_aligned)) {
    strbp = span - (aligned_read->msStart - aligned_read->read_start);
    // Check if End extends beyond the read
    if (aligned_read->msEnd - aligned_read->read_start >= span) {
      if (debug_adjust) {
	cout << "Partial alignment from left, returning false " << (strbp+diff_from_ref)- ms_length << "bp" << endl;
      }
      // Check if reasonable
      if (strbp+diff_from_ref < 0 || 
	  strbp+diff_from_ref >= (int)aligned_read->nucleotides.length() ||
	  strbp < (int)MIN_STR_LENGTH) return false;
      // Update read data
      aligned_read->partial = true;
      aligned_read->diffFromRef = (strbp+diff_from_ref)- ms_length ;
      aligned_read->detected_ms_nuc = aligned_read->reverse ? 
	aligned_read->nucleotides.substr(0,strbp+diff_from_ref) : 
	reverseComplement(aligned_read->nucleotides).substr(0,strbp+diff_from_ref);	
      return true;
    } else {
      aligned_read->partial = false;
    }
  }
  // Case 2: anchored at the right
  if ((!aligned_read->reverse & right_aligned) | (aligned_read->reverse & left_aligned)) {
    strbp = span - (start_pos+reglen-aligned_read->msEnd);
    // Check if beginning extends beyond the msStart
    if ((start_pos  + reglen - aligned_read->msStart) >= span) {
      if (debug_adjust) {
	cout << "span " << span << " start " << start_pos << " " << " end " << start_pos+reglen << " msend " << aligned_read->msEnd
	     << " msstart " << aligned_read->msStart << endl;
	cout << "Partial alignment from right, returning false " << (strbp +diff_from_ref)-ms_length<< "bp" << endl;
	cout << "strbp " << strbp << endl;
	cout << "diff " << diff_from_ref << endl;
	cout << "ms len " << ms_length << endl;
      }
      if (strbp+diff_from_ref < 0 || 
	  strbp+diff_from_ref >= (int)aligned_read->nucleotides.length() ||
	  strbp <= (int)MIN_STR_LENGTH) {
	return false;
      }
      aligned_read->partial = true;
      aligned_read->diffFromRef = (strbp +diff_from_ref)-ms_length;
      aligned_read->detected_ms_nuc = aligned_read->reverse ? 
	aligned_read->nucleotides.substr(aligned_read->nucleotides.length() - strbp-diff_from_ref,
					 strbp+diff_from_ref) :
	reverseComplement(aligned_read->nucleotides).substr(aligned_read->nucleotides.length() - 
							    strbp-diff_from_ref,strbp+diff_from_ref);	
      if (debug_adjust) {
	cerr << aligned_read->detected_ms_nuc << endl;
      }
      return true;
    } else {
      aligned_read->partial = false;
    }
  }
  return true;
}

bool BWAReadAligner::GetSTRAllele(MSReadRecord* aligned_read, const CIGAR_LIST& cigar_list, bool *partial) {
  // index where STR starts in the read
  size_t str_index = aligned_read->msStart-aligned_read->read_start +1;
  // Length of the total STR region
  size_t ms_length = aligned_read->msEnd - aligned_read->msStart -1;
  // If right end anchored
  if (*partial & (str_index < 0)) {
    str_index = 0;
    ms_length = aligned_read->msEnd - aligned_read->read_start;
  }
  // If left end anchored
  if (*partial & (str_index+ms_length >= aligned_read->nucleotides.length())) {
    ms_length = aligned_read->nucleotides.length() - str_index - 1;
  }
  // If alignment is too messy, get rid of it
  if (cigar_list.cigars.size() > MAX_CIGAR_SIZE) {
    return false;
  }
  if (cigar_list.cigars.size() == 1) {  // same as reference
    if (aligned_read->reverse) {
      aligned_read->detected_ms_nuc =
	reverseComplement(aligned_read->nucleotides).substr(str_index, ms_length);
    }
    aligned_read->detected_ms_nuc =
      aligned_read->nucleotides.substr(str_index, ms_length);
    aligned_read->diffFromRef = 0;
    return (aligned_read->detected_ms_nuc.length() > MIN_STR_LENGTH);
  }

  // get rid of the cigar score covering the left flanking region
  int str_start_in_cigar = aligned_read->msStart - aligned_read->read_start + 1;
  int pos = 0;
  int bp = 0;
  size_t cigar_index = 0;
  int diff = 0;
  int total_cigar_pos = 0;
  int diff_from_ref = 0;
  char cigar_type;
  CIGAR_LIST new_cigar_list;
  CIGAR_LIST str_cigar_list;
 
  while (pos <= str_start_in_cigar) {
    if (cigar_index >= cigar_list.cigars.size()) break;
    bp = cigar_list.cigars.at(cigar_index).num;
    cigar_type = cigar_list.cigars.at(cigar_index).cigar_type;
    if (cigar_type == 'M' || cigar_type == 'D') pos += bp;
    diff = pos - str_start_in_cigar + 1;
    if (diff >= 0) {
      if (cigar_index == 0) {
	new_cigar_list.cigars.resize(cigar_list.cigars.size() - cigar_index + 1);
	copy(cigar_list.cigars.begin()+cigar_index, cigar_list.cigars.end(), new_cigar_list.cigars.begin());
	diff -= cigar_list.cigars.at(cigar_index).num;
      } else {
	new_cigar_list.cigars.resize(cigar_list.cigars.size() - cigar_index + 2);
	copy(cigar_list.cigars.begin()+cigar_index-1, cigar_list.cigars.end(), new_cigar_list.cigars.begin());
	diff -= cigar_list.cigars.at(cigar_index-1).num;
      }
      break;
    }
    cigar_index += 1;
  }
  str_cigar_list.cigars = new_cigar_list.cigars;
  str_cigar_list.ResetString();
  new_cigar_list.cigars.clear();
  // get rid of cigar scores covering msend:end
  total_cigar_pos = ms_length;
  cigar_index = 0;
  pos = diff; // 0
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
      copy(str_cigar_list.cigars.begin(), str_cigar_list.cigars.begin()+cigar_index+1, new_cigar_list.cigars.begin());
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

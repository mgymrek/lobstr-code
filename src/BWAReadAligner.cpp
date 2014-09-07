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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "src/AlignmentUtils.h"
#include "src/BWAReadAligner.h"
#include "src/bwaseqio.h"
#include "src/nw.h"
#include "src/runtime_parameters.h"

using namespace std;

// extend reference this much to perform local realignment
const int REFEXTEND = 50;
extern unsigned char nst_nt4_table[256];
// Number of N's used to pad each reference
const int PAD = 50;
// max size of cigar score to allow
// Maximum difference between mate alignment
// and STR read alignment
const int MAX_PAIRED_DIFF = 1000;
// bwaq to use for mate trimming
const int MATE_TRIM_QUAL = 30;
// trim mate to this length
const size_t MAX_MATE_LENGTH = 100;
// maximum matches per flank
const int MAX_MAP_PER_FLANK = 1000;
// maximum mismatches for the final alignment
const int MAX_ALIGNMENT_MISMATCHES = 3;

// ** copied from BWA ** //
int64_t pos_end_multi(const bwt_multi1_t *p, int len) {
  if (p->cigar) {
    int j;
    int64_t x = p->pos;
    for (j = 0; j != p->n_cigar; ++j) {
      int op = __cigar_op(p->cigar[j]);
      if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
    }
    return x;
  } else {
    return p->pos + len;
  }
}

BWAReadAligner::BWAReadAligner(BWT* bwt_reference,
                               BNT* bnt_annotation,
                               map<int, REFSEQ>* ref_sequences,
                               gap_opt_t *opts) {
  bwase_initialize();
  _bwt_reference = bwt_reference;
  _bnt_annotation = bnt_annotation;
  _ref_sequences = ref_sequences;
  _opts = opts;
  _default_opts = gap_init_opt();
  _default_opts->max_diff = 10;
  _default_opts->max_gapo = 1;
  _default_opts->max_gape = 1;
  _default_opts->fnr = -1;
}

bool BWAReadAligner::ProcessReadPair(ReadPair* read_pair, string* err, string* messages) {
  // Initialize status variables
  read_pair->ResetAlignmentFlags();

  if (read_pair->reads.at(0).paired) {
    return ProcessPairedEndRead(read_pair, err, messages);
  } else {
    return ProcessSingleEndRead(read_pair, err, messages);
  }
}

bool BWAReadAligner::ProcessPairedEndRead(ReadPair* read_pair, string* err, string*messages) {
  // Keep XA string for multi-mappers
  string alternate_mappings = "";

  // all valid alignments for individual reads in the pair
  vector<ALIGNMENT> good_left_alignments_read1, good_right_alignments_read1,
    good_left_alignments_read2, good_right_alignments_read2;

  /* --- Step 1: Align each read separately --- */
  read_pair->read1_passed_alignment =
    ProcessRead(&read_pair->reads.at(0), read_pair->read1_passed_detection,
		&good_left_alignments_read1, &good_right_alignments_read1, err, messages);
  read_pair->read2_passed_alignment =
    ProcessRead(&read_pair->reads.at(1), read_pair->read2_passed_detection,
		&good_left_alignments_read2, &good_right_alignments_read2, err, messages);
  if (!(read_pair->read1_passed_alignment ||
	read_pair->read2_passed_alignment)) {
    return false;
  }
  if (read_pair->read1_passed_alignment) {
    read_pair->aligned_read_num = 0;
  } else {
    read_pair->aligned_read_num = 1;
  }

  /* --- Step 2: Determine if unique valid alignment and check mate pair --- */
  // Trim and align mate
  TrimMate(read_pair);
  vector<ALIGNMENT> mate_alignments;  
  if (!AlignMate(*read_pair, &mate_alignments)) {
    if (align_debug) {
      PrintMessageDieOnError("[ProcessPairedEndRead]: Didn't align mate", DEBUG);
    }
    return false;
  }
  // Set alignment info for the read that mapped
  read_pair->aligned_read_num = read_pair->read1_passed_alignment ? 0 : 1;
  const vector<ALIGNMENT>& good_left =
    read_pair->aligned_read_num == 0 ?
    good_left_alignments_read1 : good_left_alignments_read2;
  const vector<ALIGNMENT>& good_right =
    read_pair->aligned_read_num == 0 ?
    good_right_alignments_read1 : good_right_alignments_read2;
  // Find compatible alignment with mate
  ALIGNMENT matealign, final_left_alignment, final_right_alignment;
  for (size_t i = 0; i < good_left.size(); i++) {
    if (CheckMateAlignment(mate_alignments, good_left.at(i), good_right.at(i), &matealign)) {
      if (!read_pair->found_unique_alignment) {
	read_pair->found_unique_alignment = true;
	final_left_alignment = good_left.at(i);
	final_right_alignment = good_right.at(0);
      } else { // multi-mapper
	if (allow_multi_mappers) {
	  stringstream xa;
	  xa << good_left.at(i).chrom << ":" << good_left.at(i).start <<";";
	  alternate_mappings = alternate_mappings + xa.str();
	} else {
	  if (align_debug) {
	    PrintMessageDieOnError("[ProcessPairedEndRead]: Multi-mapper", DEBUG);
	  }
	  return false;
	}
      }
    }
  }
  if (!read_pair->found_unique_alignment) return false;

  /* --- Step 3: Adjust alignment and output --- */
  // try stitching first
  bool treat_as_paired = !(AlignmentUtils::StitchReads(read_pair, &final_left_alignment,
						       &final_right_alignment));
  if (OutputAlignment(read_pair, final_left_alignment, final_right_alignment,
		      matealign, alternate_mappings, treat_as_paired)) {
    return true;
  }
  return false;
}

bool BWAReadAligner::ProcessSingleEndRead(ReadPair* read_pair, string* err, string* messages) {
  // Keep track of multiple mappers
  string alternate_mappings = "";

  // all valid alignments for each read
  vector<ALIGNMENT> good_left_alignments_read1;
  vector<ALIGNMENT> good_right_alignments_read1;
  if (!ProcessRead(&read_pair->reads.at(0), read_pair->read1_passed_detection,
		   &good_left_alignments_read1,
		   &good_right_alignments_read1, err, messages)) {
    return false;
  }
  // Get rid of multi mappers
  if (good_left_alignments_read1.size() > 1) {
    if (allow_multi_mappers) {
      for (size_t i = 1; i < good_left_alignments_read1.size(); i++) {
	stringstream xa;
	xa << good_left_alignments_read1.at(i).chrom << ":" << good_left_alignments_read1.at(i).start << ";";
	alternate_mappings = alternate_mappings + xa.str();
      }
    } else {
      *err += "multi-mapper";
      return false;
    }
  }
  // Output the alignment
  read_pair->read1_passed_alignment = true;
  read_pair->found_unique_alignment = true;
  read_pair->aligned_read_num = 0;
  ALIGNMENT dummy_matealign;
  ALIGNMENT final_left_alignment = good_left_alignments_read1.front();
  ALIGNMENT final_right_alignment = good_right_alignments_read1.front();
  if (OutputAlignment(read_pair, final_left_alignment, final_right_alignment,
		      dummy_matealign, alternate_mappings, false)) {
    return true;
  } else {
    return false;
  }
  return false;
}

void BWAReadAligner::TrimMate(ReadPair* read_pair) {
  // make sure mate isn't ginormous
  if (read_pair->reads.at(1-read_pair->aligned_read_num).orig_nucleotides.size()
      > MAX_MATE_LENGTH) {
    read_pair->reads.at(1-read_pair->aligned_read_num).orig_nucleotides =
      read_pair->reads.at(1-read_pair->aligned_read_num).orig_nucleotides.substr(0, MAX_MATE_LENGTH);
    read_pair->reads.at(1-read_pair->aligned_read_num).orig_qual =
      read_pair->reads.at(1-read_pair->aligned_read_num).orig_qual.substr(0, MAX_MATE_LENGTH);
  }
  // Trim more harshly here
  string trim_nucs;
  string trim_quals;
  TrimRead(read_pair->reads.at(1-read_pair->aligned_read_num).orig_nucleotides,
	   read_pair->reads.at(1-read_pair->aligned_read_num).orig_qual,
	   &trim_nucs, &trim_quals, MATE_TRIM_QUAL);
  read_pair->reads.at(1-read_pair->aligned_read_num).orig_nucleotides = trim_nucs;
  read_pair->reads.at(1-read_pair->aligned_read_num).orig_qual = trim_quals;
}

bool BWAReadAligner::ProcessRead(MSReadRecord* read, bool passed_detection,
                                 vector<ALIGNMENT>* good_left_alignments,
                                 vector<ALIGNMENT>* good_right_alignments,
                                 string* err,
                                 string* messages) {
  if (align_debug) {
    stringstream msg;
    msg << "[ProcessRead]: ####### Processing " << read->ID << " ########";
    PrintMessageDieOnError(msg.str(), DEBUG);
  }
  if (!passed_detection) return false;

  *err = "Alignment-errors-here:";
  *messages = "Alignment-notes-here:";

  // Align the flanking regions
  bwa_seq_t* seqs = BWAAlignFlanks(*read);
  bwa_seq_t* seq_left = &seqs[0];
  bwa_seq_t* seq_right = &seqs[1];

  // fill in alignment coordinates
  vector<ALIGNMENT> left_alignments, right_alignments;
  bool one_flank_aligned = false;
  if (GetAlignmentCoordinates(seq_left, &left_alignments)) {
    one_flank_aligned++;
  }
  if (GetAlignmentCoordinates(seq_right, &right_alignments)) {
    one_flank_aligned++;
  }
  // don't need seqs anymore
  bwa_free_read_seq(2, seqs);

  // If didn't find anything, quit
  if (!one_flank_aligned) {
    *err += "No-flank-aligned;";
    return false;
  }

  // get shared alignments
  if (!GetSharedAlns(left_alignments, right_alignments,
                     good_left_alignments, good_right_alignments,
		     static_cast<int>(read->detected_ms_region_nuc.size()),
		     static_cast<int>(read->left_flank_nuc.size()),
		     static_cast<int>(read->right_flank_nuc.size()))) {
    *err += "No-shared-alignment-found;";
    if (align_debug) {
      PrintMessageDieOnError("[ProcessRead]: No shared alignment found", DEBUG);
    }
    return false;
  }

  // Set STR coordinates of shared alignments
  if (!SetSTRCoordinates(good_left_alignments, good_right_alignments, read->nucleotides.size())) {
    *err += "Could-not-set-str-coordinates;";
    return false;
  }

  // Return true if at least one good alignment
  return (good_left_alignments->size() >= 1 &&
          good_right_alignments->size() >= 1);
}

void BWAReadAligner::SetSeq(bwa_seq_t* seq, const string& flank_nuc, const string& flank_qual, const string& readid) {
  int is_comp = _opts->mode&BWA_MODE_COMPREAD;
  seq->bc[0] = 0;
  seq->tid = -1;
  seq->qual = 0;
  seq->full_len = seq->clip_len = seq->len = flank_nuc.length();
  seq->seq = reinterpret_cast<ubyte_t*>(calloc(seq->len, 1));
  seq->qual = reinterpret_cast<ubyte_t*>(calloc(seq->len, 1));
  for (int i = 0; i != seq->full_len; ++i) {
    seq->seq[i] = nst_nt4_table[static_cast<int>(flank_nuc.at(i))];
    if (fastq) {
      seq->qual[i] = (flank_qual.at(i)+33 < 126) ? flank_qual.at(i)+33 : 126;
    } else { seq->qual[i] = 126; }
  }
  seq->rseq = reinterpret_cast<ubyte_t*>(calloc(seq->full_len, 1));
  memcpy(seq->rseq, seq->seq, seq->len);
  seq_reverse(seq->len, seq->seq, 0);
  seq_reverse(seq->len, seq->rseq, is_comp);
  seq_reverse(seq->len, seq->qual, 0);
  seq->name = strdup((const char*)readid.c_str());
}

bwa_seq_t* BWAReadAligner::BWAAlignFlanks(const MSReadRecord& read) {
  // set up bwa
  bwa_seq_t *seqs, *seq_left, *seq_right;
  seqs = reinterpret_cast<bwa_seq_t*>(calloc(2, sizeof(bwa_seq_t)));
  seq_left = &seqs[0];
  seq_right = &seqs[1];

  // left in for historical reasons since BWA reverses sequence
  const string& left_flank_nuc = reverseComplement(read.left_flank_nuc);
  const string& right_flank_nuc = reverseComplement(read.right_flank_nuc);
  const string& left_qual =  reverse(read.quality_scores.
                                     substr(0, left_flank_nuc.length()));
  const string& right_qual = reverse(read.quality_scores.
                                     substr(read.quality_scores.size()
                                            - right_flank_nuc.length(),
                                            right_flank_nuc.length()));
  if (align_debug) {
    stringstream msg;
    msg << "[BWAAlignFlanks]: left flank " << left_flank_nuc << " right flank " << right_flank_nuc;
    PrintMessageDieOnError(msg.str(), DEBUG);
  }
  if (left_qual.size() != left_flank_nuc.size() || right_qual.size() != right_flank_nuc.size()) {
    PrintMessageDieOnError("[BWAAlignFlanks]: Internal error: Qual size does not match nuc size", WARNING);
    return seqs;
  }

  // Set flanking regions
  SetSeq(seq_left, left_flank_nuc, left_qual, read.ID);
  SetSeq(seq_right, right_flank_nuc, right_qual, read.ID);

  // call bwa with appropriate options (separate for each flank)
  if (seq_left->len >= min_length_to_allow_mismatches) {
    _opts->fnr = fpr;
    _opts->max_diff = max_mismatch;
    _opts->max_gapo = gap_open;
    _opts->max_gape = gap_extend;
  } else {
    _opts->fnr = -1;
    _opts->max_diff = 0;
    _opts->max_gapo = 0;
    _opts->max_gape = 0;
  }
  bwa_cal_sa_reg_gap(0, _bwt_reference->bwt,
                     1, seq_left, _opts);
  if (seq_right->len >= min_length_to_allow_mismatches) {
    _opts->fnr = fpr;
    _opts->max_diff = max_mismatch;
    _opts->max_gapo = gap_open;
    _opts->max_gape = gap_extend;
  } else {
    _opts->fnr = -1;
    _opts->max_diff = 0;
    _opts->max_gapo = 0;
    _opts->max_gape = 0;
  }
  bwa_cal_sa_reg_gap(0, _bwt_reference->bwt,
                     1, seq_right, _opts);
  return seqs;
}

void BWAReadAligner::ParseRefid(const string& refstring, ALIGNMENT* refid) {
  vector<string> items;
  split(refstring, '$', items);
  refid->id = atoi(items.at(0).c_str());
  refid->chrom = items.at(1);
  refid->start = atoi(items.at(2).c_str())+extend;
  refid->end = atoi(items.at(3).c_str())-extend;
}

bool BWAReadAligner::GetAlignmentCoordinates(bwa_seq_t* aligned_seqs,
                                             vector<ALIGNMENT>* alignments) {
  alignments->clear();
  // fill in alignment properties
  bwa_aln2seq_core(aligned_seqs->n_aln, aligned_seqs->aln,
                   aligned_seqs, 0, MAX_MAP_PER_FLANK); // last arg gives max num multi-mappers
  bwa_cal_pac_pos_core(0, _bwt_reference->bwt[1-aligned_seqs->strand],
                       &aligned_seqs[0], _opts->max_diff, _opts->fnr);

  // get coords
  for (int i = 0; i < aligned_seqs->n_multi; ++i) {
    bwt_multi1_t *q = aligned_seqs->multi + i;
    // set pos
    if (!q->strand) {
      q->pos = _bwt_reference->bwt[1]->seq_len -
        bwt_sa(_bwt_reference->bwt[1], q->pos) +
        aligned_seqs[0].len;
    } else {
      q->pos = bwt_sa(_bwt_reference->bwt[0], q->pos);
    }

    int j, seqid;
    j = pos_end_multi(q, aligned_seqs->len) - q->pos;
    if (q->pos >= _bnt_annotation->bns->l_pac) {
      continue;
    }
    bns_coor_pac2real(_bnt_annotation->bns, q->pos, j, &seqid);

    // fill in alignment info
    ALIGNMENT refid;
    ParseRefid(_bnt_annotation->bns->anns[seqid].name, &refid);
    refid.strand = q->strand;
    if (q->strand) {
      refid.pos = (refid.start-extend-PAD) +
        static_cast<int>(q->pos - _bnt_annotation->bns->anns[seqid].offset) + 1;
    } else {
      refid.pos = static_cast<int>(q->pos - _bnt_annotation->bns->anns[seqid].offset - 2*j) + (refid.start-extend-PAD);
    }
    refid.endpos = refid.pos + aligned_seqs->len;
    alignments->push_back(refid);
  }
  return (alignments->size() >= 1);
}

bool BWAReadAligner::GetSharedAlns(const vector<ALIGNMENT>& map1,
                                   const vector<ALIGNMENT>& map2,
                                   vector<ALIGNMENT>* left_refids,
                                   vector<ALIGNMENT>* right_refids,
				   int ms_len, int left_flank_len,
				   int right_flank_len) {
  if (align_debug) {
    stringstream msg;
    msg << "[GetSharedAln]: " << allow_one_flank_align << " map 1 size " << map1.size() << " map 2 size " << map2.size();
    PrintMessageDieOnError(msg.str(), DEBUG);
  }
  if (!allow_one_flank_align && (map1.size() == 0 || map2.size() == 0)) return false;
  if (map1.size() == 0 && map2.size() == 0) return false;

  if (map1.size() == 0) { // only right flank aligned
    if (!allow_multi_mappers && map2.size() > 1) return false;
    // For each right alignment, make a dummy left alignment
    for (vector<ALIGNMENT>::const_iterator it = map2.begin();
	 it != map2.end(); ++it) {
      ALIGNMENT dummy_lalign;
      dummy_lalign.id = it->id;
      dummy_lalign.left = !it->left;
      dummy_lalign.chrom = it->chrom;
      dummy_lalign.copynum = it->copynum;
      dummy_lalign.strand = it->strand;
      dummy_lalign.repeat = it->repeat;
      if (!it->strand) {
	dummy_lalign.pos = it->pos + right_flank_len + ms_len;
	dummy_lalign.endpos = it->pos + left_flank_len + ms_len + right_flank_len;
      } else {
	dummy_lalign.pos = it->pos - ms_len - left_flank_len; // approx
	dummy_lalign.endpos = it->pos - ms_len;
      }
      left_refids->push_back(dummy_lalign);
      right_refids->push_back(*it);
      if (align_debug) {
	stringstream msg;
	msg << "Dummy left_align, chrom:pos " << dummy_lalign.chrom << ":" << dummy_lalign.pos << "-" << dummy_lalign.endpos << " right pos " << it->pos << " left flank len " << left_flank_len << " right flank len " << right_flank_len <<" ms len " << ms_len << " strand " << it->strand << " id " << dummy_lalign.id;
	PrintMessageDieOnError(msg.str(), DEBUG);
      }
    }
    return true;
  } else if (map2.size() == 0) { // only left flank aligned
    if (!allow_multi_mappers && map1.size() > 1) return false;
    for (vector<ALIGNMENT>::const_iterator it = map1.begin();
	 it != map1.end(); ++it) {
      ALIGNMENT dummy_ralign;
      dummy_ralign.id = it->id;
      dummy_ralign.left = !it->left;
      dummy_ralign.chrom = it->chrom;
      dummy_ralign.copynum = it->copynum;
      dummy_ralign.strand = it->strand;
      dummy_ralign.repeat = it->repeat;
      if (!it->strand) {
	dummy_ralign.pos = it->pos - ms_len - left_flank_len;
	dummy_ralign.endpos = it->pos - ms_len;
      } else {
	dummy_ralign.pos = it->pos + left_flank_len + ms_len; // approx
	dummy_ralign.endpos = it->pos + left_flank_len + ms_len + right_flank_len;
      }
      left_refids->push_back(*it);
      right_refids->push_back(dummy_ralign);
      if (align_debug) {
	stringstream msg;
	msg << "[GetSharedAln]: Dummy right_align, chrom:pos " << dummy_ralign.chrom << ":" <<dummy_ralign.pos << " left pos " << it->pos << " left flank len " << left_flank_len << " ms len " << ms_len << " left strand " << it->strand << " id " << dummy_ralign.id;
	PrintMessageDieOnError(msg.str(), DEBUG);
      }
    }
    return true;
  } else { // both aligned, find matching ones
    // get keys that are shared and consistent
    map<int, ALIGNMENT> left_id_to_ref;
    map<int, ALIGNMENT> right_id_to_ref;
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
    // Iterate through shortest map first
    const map<int, ALIGNMENT>& smaller_id_to_ref = (right_id_to_ref.size() < left_id_to_ref.size()) ?
      right_id_to_ref : left_id_to_ref;
    const map<int, ALIGNMENT>& larger_id_to_ref = (right_id_to_ref.size() < left_id_to_ref.size()) ?
      left_id_to_ref : right_id_to_ref;
    for (map<int, ALIGNMENT>::const_iterator it = smaller_id_to_ref.begin();
	 it != smaller_id_to_ref.end(); it++) {
      int ref_key = it->first;
      if (align_debug) {
	stringstream msg;
	msg << "[GetSharedAln]: Checking smaller key " << ref_key;
	PrintMessageDieOnError(msg.str(), DEBUG);
      }
      if (larger_id_to_ref.find(ref_key) != larger_id_to_ref.end()) {
	if (right_id_to_ref[ref_key].strand == left_id_to_ref[ref_key].strand) {
	  left_refids->push_back(left_id_to_ref[ref_key]);
	  right_refids->push_back(right_id_to_ref[ref_key]);
	  if (align_debug) {
	    stringstream msg;
	    msg << "[GetSharedAln]: Found match";
	    PrintMessageDieOnError(msg.str(), DEBUG);
	  }
	  found = true;
	}
      }
    }
    return found;
  }
}

void BWAReadAligner::GetSpannedSTRs(const ALIGNMENT& lalign, const ALIGNMENT& ralign, const int& refid,
				    vector<ReferenceSTR>* spanned_ref_strs, vector<string>* repseq, size_t read_length) {
  // Get min and max coordinates
  int mincoord = lalign.pos < ralign.pos ? lalign.pos : ralign.pos;
  int maxcoord = lalign.endpos > ralign.endpos ? lalign.endpos : ralign.endpos;
  if (align_debug) {
    stringstream msg;
    msg << "[GetSpannedSTRs]: Mincoord: " << mincoord << " Maxcoord: " << maxcoord;
    PrintMessageDieOnError(msg.str(), DEBUG);
  }
  if (maxcoord - mincoord > static_cast<int>(read_length) + max_diff_ref) {
    if (align_debug) {
      stringstream msg;
      msg << "[GetSpannedSTRs]: span too large. quitting";
      PrintMessageDieOnError(msg.str(), DEBUG);
    }
    return;
  }
  const map<string, vector<ReferenceSTR> > ref_strs = (*_ref_sequences)[refid].ref_strs;
  for (map<string, vector<ReferenceSTR> >::const_iterator it = ref_strs.begin(); it != ref_strs.end(); it++) {
    const string& motif = it->first;
    for (size_t i = 0; i < it->second.size(); i++) {
      const ReferenceSTR& ref = it->second.at(i);
      if (mincoord < ref.start && maxcoord > ref.stop) {
	spanned_ref_strs->push_back(ref);
	repseq->push_back(motif);
      }
    }
  }
}

bool BWAReadAligner::SetSTRCoordinates(vector<ALIGNMENT>* good_left_alignments,
				       vector<ALIGNMENT>* good_right_alignments,
				       size_t read_length) {
  if (align_debug) {
    PrintMessageDieOnError("[SetSTRCoordinates]: Setting STR coordinates", DEBUG);
  }
  vector<ALIGNMENT> processed_left_alignments;
  vector<ALIGNMENT> processed_right_alignments;
  for (size_t i = 0; i < good_left_alignments->size(); i++) {
    const ALIGNMENT& lalign = good_left_alignments->at(i);
    const ALIGNMENT& ralign = good_right_alignments->at(i);
    int refid = lalign.id;
    // Which STR is spanned by these two alignments
    vector<ReferenceSTR> spanned_ref_strs;
    vector<string> repseqs;
    GetSpannedSTRs(lalign, ralign, refid, &spanned_ref_strs, &repseqs, read_length);
    if (align_debug) {
      stringstream msg;
      msg << "[SetSTRCoordinates]: Found " << spanned_ref_strs.size() << " spanned STRs";
      PrintMessageDieOnError(msg.str(), DEBUG);
    }
    if (spanned_ref_strs.size() == 0) continue;
    ALIGNMENT l_spanning_alignment, r_spanning_alignment;
    for (size_t j = 0; j < spanned_ref_strs.size(); j++) {
      if (j >= 1) {
	// If a single read spans multiple STRs, keep track of the other ones
	stringstream other_spanned_strs;
	other_spanned_strs << spanned_ref_strs.at(j).chrom << ":"
			   << spanned_ref_strs.at(j).start << ":"
			   << repseqs.at(j) << ";";
	l_spanning_alignment.other_spanned_strs += other_spanned_strs.str();
	r_spanning_alignment.other_spanned_strs += other_spanned_strs.str();
	continue;
      }
      float copynum = static_cast<float>(spanned_ref_strs.at(j).stop-spanned_ref_strs.at(j).start+1)/
	static_cast<float>(repseqs.at(j).size());
      copynum = ceilf(copynum * 10) / 10; // Round to 1 digit to match reference
      l_spanning_alignment.id = refid;
      l_spanning_alignment.chrom = spanned_ref_strs.at(j).chrom;
      l_spanning_alignment.start = spanned_ref_strs.at(j).start;
      l_spanning_alignment.end = spanned_ref_strs.at(j).stop;
      l_spanning_alignment.repeat = repseqs.at(j);
      l_spanning_alignment.strand = lalign.strand;
      l_spanning_alignment.pos = lalign.pos;
      l_spanning_alignment.copynum = copynum;
      r_spanning_alignment.id = refid;
      r_spanning_alignment.chrom = spanned_ref_strs.at(j).chrom;
      r_spanning_alignment.start = spanned_ref_strs.at(j).start;
      r_spanning_alignment.end = spanned_ref_strs.at(j).stop;
      r_spanning_alignment.repeat = repseqs.at(j);
      r_spanning_alignment.strand = ralign.strand;
      r_spanning_alignment.pos = ralign.pos;
      r_spanning_alignment.copynum = copynum;
      if (l_spanning_alignment.pos < r_spanning_alignment.pos) {
	l_spanning_alignment.left = true;
	r_spanning_alignment.left = false;
      } else {
	l_spanning_alignment.left = false;
	r_spanning_alignment.left = true;
      }
    }
    processed_left_alignments.push_back(l_spanning_alignment);
    processed_right_alignments.push_back(r_spanning_alignment);
  }
  *good_left_alignments = processed_left_alignments;
  *good_right_alignments = processed_right_alignments;
  return good_left_alignments->size() >= 1;
}

bool BWAReadAligner::AlignMate(const ReadPair& read_pair,
                               vector<ALIGNMENT>* mate_alignments) {
  const int& num_aligned_read = read_pair.aligned_read_num;
  const string& nucs = reverseComplement(read_pair.reads.
                                         at(1-num_aligned_read).
                                         orig_nucleotides);
  const string& qual = reverse(read_pair.reads.
                               at(1-num_aligned_read).orig_qual);

  // set up BWA alignment
  bwa_seq_t *seq = reinterpret_cast<bwa_seq_t*>(calloc(1, sizeof(bwa_seq_t)));
  SetSeq(seq, nucs, qual, read_pair.reads.at(1-num_aligned_read).ID.c_str());

  // call bwa with appropriate options
  bwa_cal_sa_reg_gap(0, _bwt_reference->bwt, 1, seq, _default_opts);

  if (seq->n_aln == 0) {
    bwa_free_read_seq(1, seq);
    return false;
  }

  // Check alignment coordinates
  if (!GetAlignmentCoordinates(seq, mate_alignments)) {
    bwa_free_read_seq(1, seq);
    return false;
  }
  bwa_free_read_seq(1, seq);
  return true;
}

bool BWAReadAligner::CheckMateAlignment(const vector<ALIGNMENT>&
                                        mate_alignments,
                                        const ALIGNMENT& left_alignment,
                                        const ALIGNMENT& right_alignment,
                                        ALIGNMENT* mate_alignment) {
  // For each, check against STR alignment
  for (vector<ALIGNMENT>::const_iterator it = mate_alignments.begin();
       it != mate_alignments.end(); ++it) {
    // find start of STR read
    const int& str_pos = (left_alignment.pos < right_alignment.pos) ?
      left_alignment.pos : right_alignment.pos;
    if ((it->chrom == left_alignment.chrom) &&
        (abs(it->pos-str_pos) <= MAX_PAIRED_DIFF) &&
        it->strand != left_alignment.strand) {
      *mate_alignment =  (*it);
      return true;
    }
  }
  return false;
}

bool BWAReadAligner::OutputAlignment(ReadPair* read_pair,
                                     const ALIGNMENT& left_alignment,
                                     const ALIGNMENT& right_alignment,
                                     const ALIGNMENT& mate_alignment,
				     const std::string& alternate_mappings,
                                     bool treat_as_paired) {
  const int& aligned_read_num = read_pair->aligned_read_num;
  read_pair->treat_as_paired = treat_as_paired;

  // Set info for aligned read
  read_pair->reads.at(aligned_read_num).chrom = left_alignment.chrom;
  read_pair->reads.at(aligned_read_num).strid = left_alignment.id;
  if (align_debug) {
    stringstream msg;
    msg << "Left aln start: " << left_alignment.pos << " Right aln start: " << right_alignment.pos << " ids " << left_alignment.id << " " << right_alignment.id;
    PrintMessageDieOnError(msg.str(), DEBUG);
  }
  if (left_alignment.left) {
    read_pair->reads.at(aligned_read_num).msStart = left_alignment.start;
    read_pair->reads.at(aligned_read_num).msEnd = left_alignment.end;
  } else {
    read_pair->reads.at(aligned_read_num).msStart = right_alignment.start;
    read_pair->reads.at(aligned_read_num).msEnd = right_alignment.end;
  }
  read_pair->reads.at(aligned_read_num).refCopyNum = left_alignment.copynum;
  read_pair->reads.at(aligned_read_num).reverse = !left_alignment.left;
  read_pair->reads.at(aligned_read_num).repseq = left_alignment.repeat;
  read_pair->reads.at(aligned_read_num).lStart = left_alignment.pos;
  read_pair->reads.at(aligned_read_num).lEnd = left_alignment.pos +
    read_pair->reads.at(aligned_read_num).left_flank_nuc.length();
  read_pair->reads.at(aligned_read_num).rStart = right_alignment.pos;
  read_pair->reads.at(aligned_read_num).rEnd = right_alignment.pos +
    read_pair->reads.at(aligned_read_num).right_flank_nuc.length();

  read_pair->alternate_mappings = alternate_mappings;
  read_pair->other_spanned_strs = left_alignment.other_spanned_strs;

  // checks to make sure the alignment is reasonable
  // coords make sense on forward strand
  if (!read_pair->reads.at(aligned_read_num).reverse &&
       ((read_pair->reads.at(aligned_read_num).lStart >=
          read_pair->reads.at(aligned_read_num).msStart)
         || (read_pair->reads.at(aligned_read_num).rEnd <=
            read_pair->reads.at(aligned_read_num).msEnd))) {
    return false;
  }
  // coords make sense on reverse strand
  if (read_pair->reads.at(aligned_read_num).reverse &&
      ((read_pair->reads.at(aligned_read_num).rStart >=
          read_pair->reads.at(aligned_read_num).msStart)
         || (read_pair->reads.at(aligned_read_num).lEnd <=
      read_pair->reads.at(aligned_read_num).msEnd))) {
    return false;
  }

  // coords of flanks are positive
  if ((read_pair->reads.at(aligned_read_num).rStart <0) ||
      (read_pair->reads.at(aligned_read_num).lStart <0) ||
      (read_pair->reads.at(aligned_read_num).rEnd <0) ||
      (read_pair->reads.at(aligned_read_num).lEnd < 0) ||
      (read_pair->reads.at(aligned_read_num).msStart < 0)) {
    return false;
  }

  // Adjust alignment and STR call
  if (!AdjustAlignment(&read_pair->reads.at(aligned_read_num))) {
    return false;
  }
  // Make sure alignment meets requirements
  if (unit) {
    if (read_pair->reads.at(aligned_read_num).diffFromRef %
        read_pair->reads.at(aligned_read_num).repseq.length() != 0) {
      return false;
    }
  }
  if ((abs(read_pair->reads.at(aligned_read_num).diffFromRef) >
        max_diff_ref) ||
      (read_pair->reads.at(aligned_read_num).mapq >= max_mapq)) {
    return false;
  }
  if (treat_as_paired) {
    // need to reset nucs/quals in case we trimmed them
    read_pair->reads.at(1-aligned_read_num).nucleotides =
      read_pair->reads.at(1-aligned_read_num).orig_nucleotides;
    read_pair->reads.at(1-aligned_read_num).quality_scores =
      read_pair->reads.at(1-aligned_read_num).orig_qual;
    read_pair->reads.at(1-aligned_read_num).read_start = mate_alignment.pos;
    read_pair->reads.at(1-aligned_read_num).reverse =
      !read_pair->reads.at(aligned_read_num).reverse;
    // get cigar
    CIGAR_LIST cigar_list;
    string aln_seq, ref_seq;
    int score;
    const size_t& reglen = read_pair->reads.
      at(1-aligned_read_num).nucleotides.length();
    const REFSEQ& refseq = _ref_sequences->
      at(read_pair->reads.at(aligned_read_num).strid);
    const size_t& start_pos = read_pair->reads.
      at(1-aligned_read_num).reverse ?
      mate_alignment.pos : mate_alignment.pos-1;
    string rseq;
    try {
      rseq = refseq.sequence.
	substr(start_pos - refseq.start + PAD, reglen);
    } catch(std::out_of_range & exception) {
      return false;      
    }
    const string& aseq = read_pair->reads.at(1-aligned_read_num).reverse ?
      reverseComplement(read_pair->reads.at(1-aligned_read_num).nucleotides) :
      read_pair->reads.at(1-aligned_read_num).nucleotides;
    const string& aligned_seq_quals =
      read_pair->reads.at(1-aligned_read_num).reverse ?
      reverse(read_pair->reads.at(1-aligned_read_num).quality_scores) :
      read_pair->reads.at(1-aligned_read_num).quality_scores;
    nw(aseq, rseq, aln_seq, ref_seq, false, &score, &cigar_list);
    // update qualities. For read pairs qual is sum of the two ends' mapq
    int edit;
    int mate_mapq = AlignmentUtils::GetMapq(aln_seq, ref_seq,
					    aligned_seq_quals, &edit);
    if (edit >= MAX_ALIGNMENT_MISMATCHES) {
      return false;
    }
    const int& read_mapq = read_pair->reads.at(aligned_read_num).mapq;
    read_pair->reads.at(1-aligned_read_num).mapq = mate_mapq+read_mapq;
    read_pair->reads.at(1-aligned_read_num).edit_dist = edit;
    read_pair->reads.at(aligned_read_num).mapq = mate_mapq+read_mapq;
    
    // need this to make it work out
    if (!read_pair->reads.at(1-aligned_read_num).reverse) {
      read_pair->reads.at(1-aligned_read_num).read_start--;
    }
    
    // make sure CIGAR is valid
    bool added_s;
    bool cigar_had_s;
    GenerateCorrectCigar(&cigar_list,read_pair->reads.at(1-aligned_read_num).
			 nucleotides, &added_s, &cigar_had_s);
    read_pair->reads.at(1-aligned_read_num).cigar_string =
      cigar_list.cigar_string;
    read_pair->reads.at(1-aligned_read_num).cigar =
      cigar_list.cigars;
  }
  return true;
}

bool BWAReadAligner::AdjustAlignment(MSReadRecord* aligned_read) {
  // get reference sequence
  const size_t& reglen = !aligned_read->reverse ?
    (aligned_read->rEnd - aligned_read->lStart) :
    (aligned_read->lEnd - aligned_read->rStart);
  if ((_ref_sequences->find(aligned_read->strid) ==
       _ref_sequences->end())) {
    PrintMessageDieOnError("[AdjustAlignment]: Ref id out of range. Problem with lobSTR index", ERROR);
  }
  const REFSEQ& refseq = _ref_sequences->at(aligned_read->strid);
  size_t start_pos = !aligned_read->reverse ?
    aligned_read->lStart-1 : aligned_read->rStart;
  string rseq;
  try {
    rseq = refseq.sequence.substr(start_pos - refseq.start-REFEXTEND/2 + PAD, reglen+REFEXTEND);
  } catch(std::out_of_range & exception) { 
    return false;
  }
  const string& aligned_seq = !aligned_read->reverse ?
    aligned_read->nucleotides :
    reverseComplement(aligned_read->nucleotides);
  const string& aligned_seq_quals = !aligned_read->reverse ?
    aligned_read->quality_scores :
    reverse(aligned_read->quality_scores);

  // update coords
  aligned_read->read_start = start_pos - REFEXTEND/2;
  aligned_read->read_end = start_pos + reglen;

  // Global alignment read vs. region aligned to
  string aligned_seq_sw, ref_seq_sw;
  int sw_score;
  CIGAR_LIST cigar_list;
  nw(aligned_seq, rseq, aligned_seq_sw, ref_seq_sw,
     false, &sw_score, &cigar_list);
  if (align_debug) {
    stringstream msg;
    msg << "Aseq " << aligned_seq_sw;
    PrintMessageDieOnError(msg.str(), DEBUG);
    stringstream msg2;
    msg2 << "Rseq " << ref_seq_sw;
    PrintMessageDieOnError(msg2.str(), DEBUG);
  }
  // get rid of end gaps
  if (cigar_list.cigars.at(0).cigar_type == 'D') {
    const int& num = cigar_list.cigars.at(0).num;
    aligned_read->read_start += num;
    cigar_list.cigars.erase(cigar_list.cigars.begin());
  }
  if (cigar_list.cigars.at(cigar_list.cigars.size() - 1).cigar_type == 'D') {
    cigar_list.cigars.erase(cigar_list.cigars.end() - 1);
  }
  if (cigar_list.cigars.at(0).cigar_type == 'I') {
    cigar_list.cigars.at(0).cigar_type = 'S';
  }
  if (cigar_list.cigars.at(cigar_list.cigars.size() - 1).cigar_type == 'D') {
    cigar_list.cigars.at(cigar_list.cigars.size() - 1).cigar_type = 'S';
  }
  cigar_list.ResetString();

  // Update info in aligned read
  int edit;
  aligned_read->mapq = AlignmentUtils::GetMapq(aligned_seq_sw, ref_seq_sw,
					       aligned_seq_quals, &edit);
  aligned_read->edit_dist = edit;

  // make sure CIGAR is valid
  bool added_s;
  bool cigar_had_s;
  GenerateCorrectCigar(&cigar_list, aligned_read->nucleotides, &added_s, &cigar_had_s);
  aligned_read->cigar = cigar_list.cigars;
  aligned_read->cigar_string = cigar_list.cigar_string;

  // Check if alignment is reasonably good
  if (sw_score < min_sw_score) {
    return false;
  }

  // Update STR allele
  try {
    return AlignmentUtils::GetSTRAllele(aligned_read, cigar_list);
  } catch(std::out_of_range & exception) {
    return false;
  }
  return false;
}

BWAReadAligner::~BWAReadAligner() {}

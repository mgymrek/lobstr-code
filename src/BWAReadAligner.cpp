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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "src/BWAReadAligner.h"
#include "src/bwaseqio.h"
#include "src/nw.h"
#include "src/runtime_parameters.h"

using namespace std;

// extend reference this much to perform local realignment
const int REFEXTEND = 10;
extern unsigned char nst_nt4_table[256];
// Number of N's used to pad each reference
const int PAD = 50;
// max size of cigar score to allow
// more than this is likely bad alignment
const size_t MAX_CIGAR_SIZE = 8;
// Maximum difference between mate alignment
// and STR read alignment
const int MAX_PAIRED_DIFF = 1000;
// BWA % mismatching for mate alignment
const float MATE_FNR = 0.1;
// BWA # mismatches allowed in mate alignment
const int MATE_MISMATCH = 10;
// BWA # gap opens allowed in mate alignment
const int MATE_GAPO = 1;
// BWA # gap extends allowed in mate alignment
const int MATE_GAPE = 10;
// bwaq to use for mate trimming
const int MATE_TRIM_QUAL = 30;
// trim mate to this length
const size_t MAX_MATE_LENGTH = 100;
// Minimum number of bp for stitch overlap
const size_t MIN_STITCH_OVERLAP = 16;
// Percent identity required to stitch
const float STITCH_REQUIRED_SCORE = 0.8;
// Allowed difference in score between returned stitch
// and next best stitch
const float STITCH_DIFF = 0.1;
// min allowed distance from STR boundary to read ends
size_t MIN_DIST_FROM_END = 8;

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
  _default_opts->max_diff = MATE_MISMATCH;
  _default_opts->max_gapo = MATE_GAPO;
  _default_opts->max_gape = MATE_GAPE;
  _default_opts->fnr = -1;//MATE_FNR;

  cigar_debug = false;
  stitch_debug = false;
}

bool BWAReadAligner::ProcessReadPair(ReadPair* read_pair, string* err, string* messages) {
  // Initialize status variables
  read_pair->read1_passed_alignment = false;
  read_pair->read2_passed_alignment = false;
  read_pair->found_unique_alignment = false;
  read_pair->aligned_read_num = -1;

  // Keep XA string for multi-mappers
  string alternate_mappings = "";

  if (read_pair->reads.at(0).paired) {
    // all valid alignments for individual reads in the pair
    vector<ALIGNMENT> good_left_alignments_read1;
    vector<ALIGNMENT> good_right_alignments_read1;
    vector<ALIGNMENT> good_left_alignments_read2;
    vector<ALIGNMENT> good_right_alignments_read2;

    // Step 1: Align each read separately
    if (read_pair->read1_passed_detection) {
      if (ProcessRead(&read_pair->reads.at(0),
                      &good_left_alignments_read1,
                      &good_right_alignments_read1, err, messages)) {
        read_pair->read1_passed_alignment = true;
      }
    }
    if (read_pair->read2_passed_detection) {
      if (ProcessRead(&read_pair->reads.at(1),
                      &good_left_alignments_read2,
                      &good_right_alignments_read2, err, messages)) {
        read_pair->read2_passed_alignment = true;
      }
    }
    // Need at least one read to pass
    if (!(read_pair->read1_passed_alignment ||
          read_pair->read2_passed_alignment)) {
      return false;
    }
    if (align_debug) {
      stringstream msg;
      msg << "[BWAReadAligner]: read1 passed "
           << read_pair->read1_passed_alignment;
      PrintMessageDieOnError(msg.str(), DEBUG);
      msg.clear();
      msg << "[BWAReadAligner]: read2 passed "
           << read_pair->read2_passed_alignment;
      PrintMessageDieOnError(msg.str(), DEBUG);
    }
    // Step 2: Determine if unique valid alignment
    if ((read_pair->read1_passed_alignment ||
         read_pair->read2_passed_alignment)) {
      if (align_debug) {
        PrintMessageDieOnError("[BWAReadAligner]: checking paired alignment", DEBUG);
      }
      ALIGNMENT matealign;
      size_t num_alignments = 0;
      size_t index_of_hit;

      // Get info for the read that aligned
      if (read_pair->read1_passed_alignment) {
	read_pair->aligned_read_num = 0;
	num_alignments = good_left_alignments_read1.size();
      } else {
        // only read 2 aligned
        read_pair->aligned_read_num = 1;
        num_alignments = good_left_alignments_read2.size();
      }

      // Get good alignments
      const vector<ALIGNMENT>& good_left =
        read_pair->aligned_read_num == 0 ?
        good_left_alignments_read1 :
        good_left_alignments_read2;
      const vector<ALIGNMENT>& good_right =
        read_pair->aligned_read_num == 0 ?
        good_right_alignments_read1 :
        good_right_alignments_read2;

      // If both aligned, find compatible
      if (!read_pair->found_unique_alignment &&
          read_pair->read1_passed_alignment &&
          read_pair->read2_passed_alignment) {
        if (align_debug) {
          stringstream msg;
          msg << "[BWAReadAligner]: checking for compatible alignment"
              << " read 1 size " << good_left_alignments_read1.size()
              << " read 2 size " << good_left_alignments_read2.size();
          PrintMessageDieOnError(msg.str(), DEBUG);
        }
        // find compatible alignment
        size_t index_of_mate;
        if (FindCompatibleAlignment(good_left_alignments_read1,
                                    good_left_alignments_read2,
                                    good_right_alignments_read1,
                                    good_right_alignments_read2,
                                    &index_of_hit, &index_of_mate,
				    &alternate_mappings)) {
          if (align_debug) {
            PrintMessageDieOnError("[BWAReadAligner]: found compatible alignments", DEBUG);
          }
          if (read_pair->aligned_read_num == 0) {
            matealign = read_pair->
              reads.at(1).reverse ?
              good_left_alignments_read2.at(index_of_mate) :
              good_right_alignments_read2.at(index_of_mate);
          } else {
            matealign = read_pair->
              reads.at(0).reverse ?
              good_left_alignments_read1.at(index_of_hit) :
              good_right_alignments_read1.at(index_of_hit);
            index_of_hit = index_of_mate;
          }
          read_pair->found_unique_alignment = true;
        }
      }

      // Still didn't find alignment, or only one aligned,
      // Check that mate is compatible
      if (!read_pair->found_unique_alignment) {
        vector<ALIGNMENT> mate_alignments;
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
        if (align_debug) {
          PrintMessageDieOnError("[BWAReadAligner]: checking for mate alignment " + trim_nucs + " " + trim_quals, DEBUG);
        }
        if (!AlignMate(*read_pair, &mate_alignments,
                       read_pair->reads.at(read_pair->aligned_read_num).repseq)) {
          if (align_debug) {
            PrintMessageDieOnError("[BWAReadAligner]: Mate not aligned.", DEBUG);
          }
        } else {
          for (size_t i = 0; i < num_alignments; i++) {
            const ALIGNMENT& lalign = good_left.at(i);
            const ALIGNMENT& ralign = good_right.at(i);
            if (CheckMateAlignment(mate_alignments, lalign, ralign,
                                   &matealign)) {
	      index_of_hit = i;
              if (!read_pair->found_unique_alignment) {
                read_pair->found_unique_alignment = true;
              } else {
                // multiple mapper, more than one good hit
		if (allow_multi_mappers) {
		  stringstream xa;
		  xa << lalign.chrom << ":" << lalign.start << ";";
		  alternate_mappings = alternate_mappings + xa.str();
		} else {
		  return false;
		}
              }
            } else {
	      if (align_debug) {
		PrintMessageDieOnError("Check mate alignment failed", DEBUG);
	      }
	    }
          }
        }
      }

      // Step 3: Adjust alignment and output
      if (align_debug) {
        PrintMessageDieOnError("[BWAReadAligner]: adjust paired alignment and output", DEBUG);
      }
      if (read_pair->found_unique_alignment) {
        ALIGNMENT final_left_alignment =
          good_left.at(index_of_hit);
        ALIGNMENT final_right_alignment =
          good_right.at(index_of_hit);
        // try stitching first
        if (align_debug) {
          PrintMessageDieOnError("[BWAReadAligner]: Try stitching", DEBUG);
        }
	if (StitchReads(read_pair, &final_left_alignment,
                        &final_right_alignment)) {
          if (OutputAlignment(read_pair, final_left_alignment,
			      final_right_alignment,
			      matealign, alternate_mappings, false)) {
	    return true;
	  }
        } else {
          if (OutputAlignment(read_pair, final_left_alignment,
			      final_right_alignment,
			      matealign, alternate_mappings, true)) {
	    return true;
	  }
	}
	// If we made it this far and didn't align, reset the
	// read so we don't try again
	if (align_debug) {
	  PrintMessageDieOnError("[BWAReadAligner]: Failed to align", DEBUG);
	}
	read_pair->reads.at(0).ms_repeat_next_best_period = 0;
	read_pair->reads.at(1).ms_repeat_next_best_period = 0;
	return false;
      }
    }
  } else {  // single end reads
    // all valid alignments for each read
    vector<ALIGNMENT> good_left_alignments_read1;
    vector<ALIGNMENT> good_right_alignments_read1;
    if (read_pair->read1_passed_detection) {
      if (ProcessRead(&read_pair->reads.at(0),
                      &good_left_alignments_read1,
                      &good_right_alignments_read1, err, messages)) {
        // Get rid of multi mappers
	if (allow_multi_mappers) {
	  if (good_left_alignments_read1.size() == 0) return false;
	  if (good_left_alignments_read1.size() > 1) {
	    for (size_t i = 1; i < good_left_alignments_read1.size(); i++) {
	      stringstream xa;
	      xa << good_left_alignments_read1.at(i).chrom << ":" << good_left_alignments_read1.at(i).start << ";";
	      alternate_mappings = alternate_mappings + xa.str();
	    }
	  }
	} else {
	  if (good_left_alignments_read1.size() > 1) {
	    return false;
	  }
	}
        if (align_debug) {
          PrintMessageDieOnError("[BWAReadAligner]: checking single end alignment.", DEBUG);
        }
        // Output the alignment
        read_pair->read1_passed_alignment = true;
        read_pair->found_unique_alignment = true;
        read_pair->aligned_read_num = 0;
        ALIGNMENT dummy_matealign;
        if (OutputAlignment(read_pair, good_left_alignments_read1.front(),
			    good_right_alignments_read1.front(),
			    dummy_matealign, alternate_mappings, false)) {
	  return true;
	} else {
	  // reset read so we don't try again
	  read_pair->reads.at(0).ms_repeat_next_best_period = 0;
	  return false;
	}
      }
    }
  }
  return false;
}

bool BWAReadAligner::ProcessRead(MSReadRecord* read,
                                 vector<ALIGNMENT>* good_left_alignments,
                                 vector<ALIGNMENT>* good_right_alignments,
                                 string* err,
                                 string* messages) {
  *err = "Alignment-errors-here:";
  *messages = "Alignment-notes-here:";
  if (align_debug) {
    PrintMessageDieOnError("[ProcessRead]: " + read->ID + " " + read->left_flank_nuc + " " + read->right_flank_nuc, DEBUG);
  }
  // do the flanks consist of all repeats?
  bool left_all_repeats = false;
  bool right_all_repeats = false;

  // Detected motif
  const string& repseq = read->repseq;

  // check that we have that repseq
  if (_bwt_references->find(repseq) == _bwt_references->end()) {
    *err += "Repseq-doesn't-exist-in-ref;";
    return false;
  }

  // Check if flanks are perfect repeats
  if (IsPerfectRepeat(read->left_flank_nuc, repseq) ||
      IsPerfectRepeat(read->left_flank_nuc, reverseComplement(repseq))) {
    read->left_perfect_repeat = true;
    left_all_repeats = true;
  }
  if (IsPerfectRepeat(read->right_flank_nuc, repseq) ||
      IsPerfectRepeat(read->right_flank_nuc, reverseComplement(repseq))) {
    read->right_perfect_repeat = true;
    right_all_repeats = true;
  }
  // don't continue if both are fully repetitive
  if (left_all_repeats && right_all_repeats) {
    *err += "Left-and-right-flanks-are-fully-repeitive;";
    return false;
  }

  if (align_debug) {
    PrintMessageDieOnError("[ProcessRead]: align flanks", DEBUG);
  }
  // Align the flanking regions
  if (read->left_flank_nuc.size() > read->quality_scores.size() ||
      read->right_flank_nuc.size() > read->quality_scores.size()) {
    PrintMessageDieOnError("[ProcessRead]: Internal error: quality score length invalid", WARNING);
  }
  bwa_seq_t* seqs = BWAAlignFlanks(*read);
  bwa_seq_t* seq_left = &seqs[0];
  bwa_seq_t* seq_right = &seqs[1];

  if (align_debug) {
    PrintMessageDieOnError("[ProcessRead]: fill in coords", DEBUG);
  }
  // fill in alignment coordinates
  vector<ALIGNMENT> left_alignments;
  vector<ALIGNMENT> right_alignments;
  if (!left_all_repeats) {
    if (!GetAlignmentCoordinates(seq_left, repseq, &left_alignments)) {
      bwa_free_read_seq(2, seqs);
      *err += "No-left-alignment;";
      return false;
    }
  }
  if (!right_all_repeats) {
    if (!GetAlignmentCoordinates(seq_right, repseq, &right_alignments)) {
      bwa_free_read_seq(2, seqs);
      *err += "No-right-alignment;";
      return false;
    }
  }

  // don't need seqs anymore
  bwa_free_read_seq(2, seqs);

  // get shared alignments
  vector<ALIGNMENT> left_refids;
  vector<ALIGNMENT> right_refids;
  if (!GetSharedAlns(left_alignments, right_alignments,
                     &left_refids, &right_refids)) {
    *err += "No-shared-alignment-found;";
    return false;
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
  if (align_debug) {
    PrintMessageDieOnError("[BWAAlignFlanks]: set up", DEBUG);
  }
  // set up bwa
  int is_comp = _opts->mode&BWA_MODE_COMPREAD;
  bwa_seq_t *seqs, *seq_left, *seq_right;
  seqs = reinterpret_cast<bwa_seq_t*>(calloc(2, sizeof(bwa_seq_t)));
  seq_left = &seqs[0];
  seq_right= &seqs[1];

  if (align_debug) {
    PrintMessageDieOnError("[BWAAlignFlanks]: reverse", DEBUG);
  }
  // left in for historical reasons since BWA reverses sequence
  const string& left_flank_nuc = reverseComplement(read.left_flank_nuc);
  const string& right_flank_nuc = reverseComplement(read.right_flank_nuc);

  if (align_debug) {
    PrintMessageDieOnError("[BWAAlignFlanks]: left flank", DEBUG);
  }
  // LEFT FLANKING REGION
  const string& left_qual =  reverse(read.quality_scores.
                                     substr(0, left_flank_nuc.length()));
  if (left_qual.size() != left_flank_nuc.size()) {
    PrintMessageDieOnError("[BWAAlignFlanks]: Qual size does not match nuc size", WARNING);
    return seqs;
  }
  seq_left->bc[0] = 0;
  seq_left->tid = -1;
  seq_left->qual = 0;
  seq_left->full_len = seq_left->clip_len = seq_left->len =
    left_flank_nuc.length();
  seq_left->seq = reinterpret_cast<ubyte_t*>(calloc(seq_left->len, 1));
  seq_left->qual = reinterpret_cast<ubyte_t*>(calloc(seq_left->len, 1));
  for (int i = 0; i != seq_left->full_len; ++i) {
    seq_left->seq[i] =
      nst_nt4_table[static_cast<int>(left_flank_nuc.at(i))];
    if (fastq) {
      seq_left->qual[i] = left_qual.at(i)+33 < 126 ?
        left_qual.at(i) + 33:126;
    } else {
      seq_left->qual[i] = 126;
    }
  }
  seq_left->rseq = reinterpret_cast<ubyte_t*>(calloc(seq_left->full_len, 1));
  memcpy(seq_left->rseq, seq_left->seq, seq_left->len);
  seq_reverse(seq_left->len, seq_left->seq, 0);
  seq_reverse(seq_left->len, seq_left->rseq, is_comp);
  seq_reverse(seq_left->len, seq_left->qual, 0);
  seq_left->name = strdup((const char*)read.ID.c_str());

  if (align_debug) {
    PrintMessageDieOnError("[BWAAlignFlanks]: right flank", DEBUG);
  }
  // RIGHT FLANKING REGION
  const string& right_qual = reverse(read.quality_scores.
                                     substr(read.quality_scores.size()
                                            - right_flank_nuc.length(),
                                            right_flank_nuc.length()));
  seq_right->bc[0] = 0;
  seq_right->tid = -1;
  seq_right->qual = 0;
  seq_right->full_len = seq_right->clip_len = seq_right->len =
    right_flank_nuc.length();
  seq_right->seq = reinterpret_cast<ubyte_t*>(calloc(seq_right->len, 1));
  seq_right->qual = reinterpret_cast<ubyte_t*>(calloc(seq_right->len, 1));
  for (int i = 0; i != seq_right->full_len; ++i) {
    seq_right->seq[i] =
      nst_nt4_table[static_cast<int>(right_flank_nuc.at(i))];
    if (fastq) {
      seq_right->qual[i] = right_qual.at(i)+33 < 126 ?
        right_qual.at(i) + 33:126;
    } else {
      seq_right->qual[i] = 126;
    }
  }
  seq_right->rseq = reinterpret_cast<ubyte_t*>(calloc(seq_right->full_len, 1));
  memcpy(seq_right->rseq, seq_right->seq, seq_right->len);
  seq_reverse(seq_right->len, seq_right->seq, 0);
  seq_reverse(seq_right->len, seq_right->rseq, is_comp);
  seq_reverse(seq_right->len, seq_right->qual, 0);
  seq_right->name = strdup((const char*)read.ID.c_str());

  // call bwa with appropriate options
  bwa_cal_sa_reg_gap(0, _bwt_references->at(read.repseq).bwt,
                     2, seqs, _opts);

  return seqs;
}

void BWAReadAligner::ParseRefid(const string& refstring, ALIGNMENT* refid) {
  vector<string> items;
  split(refstring, '$', items);
  refid->id = atoi(items.at(0).c_str());
  refid->left = (items.at(1) == "L");
  refid->chrom = items.at(2);
  refid->start = atoi(items.at(3).c_str())+extend;
  refid->end = atoi(items.at(4).c_str());
  refid->copynum = atof(items.at(6).c_str());
  if (items.size() > 7)
    refid->name = items.at(7);
}

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


bool BWAReadAligner::GetAlignmentCoordinates(bwa_seq_t* aligned_seqs,
                                             const std::string&repseq,
                                             vector<ALIGNMENT>* alignments) {
  // fill in alignment properties
  bwa_aln2seq_core(aligned_seqs->n_aln, aligned_seqs->aln,
                   aligned_seqs, 0, 1000);
  bwa_cal_pac_pos_core(0, _bwt_references->
                       at(repseq).bwt[1-aligned_seqs->strand],
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
        static_cast<int>(q->pos - _bnt_annotations->
                         at(repseq).bns->anns[seqid].offset) + 1;
    } else {
      refid.pos = static_cast<int>(q->pos - _bnt_annotations->at(repseq).
                                   bns->anns[seqid].offset - 2*j) +
        (refid.start-extend-PAD);
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
  for (map<int, ALIGNMENT>::const_iterator it = right_id_to_ref.begin();
       it != right_id_to_ref.end(); ++it) {
    int ref_key = it->first;
    if (left_id_to_ref.find(ref_key) != left_id_to_ref.end()) {
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
  const string& nucs = reverseComplement(read_pair.reads.
                                         at(1-num_aligned_read).
                                         orig_nucleotides);
  const string& qual = reverse(read_pair.reads.
                               at(1-num_aligned_read).orig_qual);

  // set up BWA alignment
  int is_comp = _opts->mode&BWA_MODE_COMPREAD;
  bwa_seq_t *seq = reinterpret_cast<bwa_seq_t*>(calloc(1, sizeof(bwa_seq_t)));
  seq->bc[0] = 0;
  seq->tid = -1;
  seq->qual = 0;
  seq->full_len = seq->clip_len = seq->len = nucs.length();
  seq->seq = reinterpret_cast<ubyte_t*>(calloc(seq->len, 1));
  seq->qual = reinterpret_cast<ubyte_t*>(calloc(seq->len, 1));
  for (int i = 0; i != seq->full_len; ++i) {
    seq->seq[i] = nst_nt4_table[static_cast<int>(nucs.at(i))];
    seq->qual[i] = qual.at(i)+33 < 126 ? qual.at(i) + 33:126;
  }
  seq->rseq = reinterpret_cast<ubyte_t*>(calloc(seq->full_len, 1));
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

bool BWAReadAligner::CheckMateAlignment(const vector<ALIGNMENT>&
                                        mate_alignments,
                                        const ALIGNMENT& left_alignment,
                                        const ALIGNMENT& right_alignment,
                                        ALIGNMENT* mate_alignment) {
  // For each, check against STR alignment
  for (vector<ALIGNMENT>::const_iterator it = mate_alignments.begin();
       it != mate_alignments.end(); ++it) {
    // Check chrom, position, and strand
    if (align_debug) {
      stringstream msg;
      msg << "[CheckMateAlignment]: check mate finds "
          << it->chrom << " "
          << (abs(it->pos-left_alignment.pos)) << " "
          << it->strand;
      PrintMessageDieOnError(msg.str(), DEBUG);
    }
    // find start of STR read
    const int& str_pos = left_alignment.left ?
      left_alignment.pos : right_alignment.pos;
    // check chrom, pos, and strand
    if (align_debug) {
      stringstream msg;
      msg << "[CheckMateAlignment]: checks "
	  << it->chrom << ": " << left_alignment.chrom
	  << " " << abs(it->pos-str_pos) << ": " << MAX_PAIRED_DIFF
	  << " " << it->strand << ": " << left_alignment.strand;
      PrintMessageDieOnError(msg.str(), DEBUG);
    }
    if ((it->chrom == left_alignment.chrom) &&
        (abs(it->pos-str_pos) <= MAX_PAIRED_DIFF) &&
        it->strand != left_alignment.strand) {
      *mate_alignment =  (*it);
      return true;
    }
  }
  return false;
}

bool BWAReadAligner::FindCompatibleAlignment(const vector<ALIGNMENT>&
                                             good_left1,
                                             const vector<ALIGNMENT>&
                                             good_left2,
                                             const vector<ALIGNMENT>&
                                             good_right1,
                                             const vector<ALIGNMENT>&
                                             good_right2,
                                             size_t* index_of_hit,
                                             size_t* index_of_mate,
					     string* alternate_mappings) {
  if (align_debug) {
    PrintMessageDieOnError("[FindCompatibleAlignment]: Looking for compatible alignment", DEBUG);
  }
  bool found_unique = false;
  // For each good left 1, check good left 2
  for (size_t i1 = 0; i1 < good_left1.size(); i1++) {
    for (size_t i2 = 0; i2 < good_left2.size(); i2++) {
      const ALIGNMENT& l1 = good_left1.at(i1);
      const ALIGNMENT& l2 = good_left2.at(i2);
      if ((abs(l1.pos-l2.pos) <= MAX_PAIRED_DIFF) &&
          (l1.strand != l2.strand) &&
          (l1.chrom == l2.chrom)) {
        if (found_unique) {
	  if (allow_multi_mappers) {
	    stringstream xa;
	    xa << l1.chrom << ":" << l1.start << ";";
	    *alternate_mappings = *alternate_mappings + xa.str();
	  } else {
	    return false;
	  }
        } else {
          found_unique = true;
          *index_of_hit = i1;
          *index_of_mate = i2;
        }
      }
    }
  }
  return found_unique;
}


bool BWAReadAligner::StitchReads(ReadPair* read_pair,
                                 ALIGNMENT* left_alignment,
                                 ALIGNMENT* right_alignment) {
  try {
    if (align_debug) {
      PrintMessageDieOnError("[StitchReads]: stitching reads", DEBUG);
    }
    // Set up
    const int& num_aligned_read = read_pair->aligned_read_num;
    string seq1 = read_pair->reads.at(num_aligned_read).orig_nucleotides;
    string seq2 = reverseComplement(read_pair->reads.
                                    at(1-num_aligned_read).
                                    orig_nucleotides);
    string seq1_qual = read_pair->reads.at(num_aligned_read).orig_qual;
    string seq2_qual = reverse(read_pair->reads.
                               at(1-num_aligned_read).orig_qual);
    bool best_stitch_is_backwards = false;
    if (stitch_debug || align_debug) {
      PrintMessageDieOnError("[StitchReads]: seq1 " + seq1 + " seq2 " + seq2, DEBUG);
    }
    vector<float> scores;
    scores.push_back(0);
    float score, max_score = 0;
    size_t overlap_len, max_score_index = -1;
    // Gradually bring ends together and try to stitch
    for (size_t i = 0; i <= seq1.length() - MIN_STITCH_OVERLAP; i++) {
      score = 0;
      overlap_len = seq1.length()-i;
      for (size_t j = 0; j < overlap_len; j++) {
        if (j >=  seq2.length()) {
          score = 0;
        } else {
          if (seq1.at(i+j) == seq2.at(j)) {
            score += 1;
          }
        }
      }
      if (score/overlap_len >= max_score) {
        max_score = score/overlap_len;
        max_score_index = i;
      }
      scores.push_back(score/overlap_len);
    }
    // Other direction
    for (size_t i = 0; i <= seq2.length() - MIN_STITCH_OVERLAP; i++) {
      score = 0;
      overlap_len = seq2.length()-i;
      for (size_t j = 0; j < overlap_len; j++) {
        if (j >=  seq1.length()) {
          score = 0;
        } else {
          if (seq2.at(i+j) == seq1.at(j)) {
            score += 1;
          }
        }
      }
      if (score/overlap_len >= max_score) {
        max_score = score/overlap_len;
        max_score_index = i + (seq1.length() - MIN_STITCH_OVERLAP) + 1;
      }
      scores.push_back(score/overlap_len);
    }
    if (stitch_debug) {
      PrintMessageDieOnError("[StitchReads]: checking for too many matches", DEBUG);
    }

    // Check if too many matches
    for (size_t i = 0; i < scores.size(); i++) {
      if ((max_score - scores.at(i) <= STITCH_DIFF) && i != max_score_index+1) {
        if (stitch_debug) {
          PrintMessageDieOnError("[StitchReads]: Returning false, too many matches i", DEBUG);
        }
        return false;
      }
    }
    if (max_score_index >= (seq1.length()-MIN_STITCH_OVERLAP)) {
      best_stitch_is_backwards = true;
      max_score_index = max_score_index - seq1.length()
        + MIN_STITCH_OVERLAP - 1;
      string tmp = seq1;
      seq1 = seq2;
      seq2 = tmp;
      tmp = seq1_qual;
      seq1_qual = seq2_qual;
      seq2_qual = tmp;
    }

    if (stitch_debug) {
      PrintMessageDieOnError("[StitchReads]: Checking that stitch is good enough", DEBUG);
    }

    // Check if stitch is good enough
    overlap_len = seq1.length() - max_score_index - 1;
    if ((overlap_len < MIN_STITCH_OVERLAP) ||
        (max_score < STITCH_REQUIRED_SCORE)) {
      return false;
    }
    string stitched_string = seq1.
      substr(0, static_cast<int>(max_score_index));
    string stitched_qual = seq1_qual.
      substr(0, static_cast<int>(max_score_index));
    string na, nb, qa, qb;

    for (size_t i = 0; i <= overlap_len; i++) {
      na = seq1.substr(max_score_index+i, 1);
      nb = seq2.substr(i, 1);
      qa = seq1_qual.substr(max_score_index+i, 1);
      qb = seq2_qual.substr(i, 1);
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
    stitched_string.append(seq2.substr(overlap_len + 1));
    stitched_qual.append(seq2_qual.substr(overlap_len + 1));

    if (stitch_debug) {
      PrintMessageDieOnError("[StitchReads]: orig string 1 " + seq1, DEBUG);
      PrintMessageDieOnError("[StitchReads]: orig string 2 " + seq2, DEBUG);
      PrintMessageDieOnError("[StitchReads]: orig qual 1 " + seq1_qual, DEBUG);
      PrintMessageDieOnError("[StitchReads]: orig qual 2 " + seq2_qual, DEBUG);
      PrintMessageDieOnError("[StitchReads]: stitched str  " + stitched_string, DEBUG);
    }
    // put stitched info in aligned read
    if (stitch_debug || align_debug) {
      PrintMessageDieOnError("[StitchReads]: " + stitched_string, DEBUG);
      PrintMessageDieOnError("[StitchReads]: " + stitched_qual, DEBUG);
    }
    read_pair->reads.at(num_aligned_read).nucleotides = stitched_string;
    read_pair->reads.at(num_aligned_read).quality_scores = stitched_qual;
    read_pair->reads.at(num_aligned_read).right_flank_nuc =
      stitched_string.substr(seq1.length() - read_pair->
                             reads.at(num_aligned_read).
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
  } catch(std::out_of_range & exception) {
    if (align_debug) {
      PrintMessageDieOnError("[StitchReads]: stitching failed. Substring out of range error.", DEBUG);
    }
    return false;
  }
}

bool BWAReadAligner::OutputAlignment(ReadPair* read_pair,
                                     const ALIGNMENT& left_alignment,
                                     const ALIGNMENT& right_alignment,
                                     const ALIGNMENT& mate_alignment,
				     const std::string& alternate_mappings,
                                     bool treat_as_paired) {
  if (align_debug) {
    PrintMessageDieOnError("[BWAReadAligner]: Output alignment", DEBUG);
  }
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

  read_pair->alternate_mappings = alternate_mappings;

  if (align_debug) {
    PrintMessageDieOnError("[BWAReadAligner]: Checkalignment", DEBUG);
  }

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
    if (align_debug) {
      PrintMessageDieOnError("[BWAReadAligner]: Discarding: reverse " \
                             "strand coords don't make sense.", DEBUG);
    }
    return false;
  }

  // coords of flanks are positive
  if ((read_pair->reads.at(aligned_read_num).rStart <0) ||
      (read_pair->reads.at(aligned_read_num).lStart <0) ||
      (read_pair->reads.at(aligned_read_num).rEnd <0) ||
      (read_pair->reads.at(aligned_read_num).lEnd < 0) ||
      (read_pair->reads.at(aligned_read_num).msStart < 0)) {
    if (align_debug) {
      PrintMessageDieOnError("[BWAReadAligner]: Discarding: negative coords found.", DEBUG);
    }
    return false;
  }

  if (align_debug) {
    PrintMessageDieOnError("[BWAReadAligner]: Adjusting alignment", DEBUG);
  }
  
  // Adjust alignment and STR call
  if (!AdjustAlignment(&read_pair->reads.at(aligned_read_num))) {
    if (align_debug) {
      PrintMessageDieOnError("[BWAReadAligner]: Returning false: "	\
			     " AdjustAlignment failed.", DEBUG);
    }
    return false;
  }
  if (align_debug) {
    PrintMessageDieOnError("[BWAReadAligner]: Checking unit requirements", DEBUG);
  }
  // Make sure alignment meets requirements
  if (unit) {
    if (read_pair->reads.at(aligned_read_num).diffFromRef %
        read_pair->reads.at(aligned_read_num).ms_repeat_best_period != 0) {
      if (align_debug) {
        PrintMessageDieOnError("[BWAReadAligner]: returning false (unit failed)", DEBUG);
      }
      return false;
    }
  }
  if (align_debug) {
    PrintMessageDieOnError("[BWAReadAligner]: Checking diff and mapq requirements", DEBUG);
  }
  if ((abs(read_pair->reads.at(aligned_read_num).diffFromRef) >
        max_diff_ref) ||
      (read_pair->reads.at(aligned_read_num).mapq >= max_mapq)) {
    if (align_debug) {
      PrintMessageDieOnError("[BWAReadAligner]: returning false " \
                             "(maxdiffref or mapq fail)", DEBUG);
    }
    return false;
  }
  if (align_debug) {
    PrintMessageDieOnError("[BWAReadAligner]: Checking treat as paired", DEBUG);
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
    if (align_debug) {
      PrintMessageDieOnError("[BWAReadAligner]: Checking get CIGAR", DEBUG);
    }
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
	substr(start_pos - refseq.start, reglen);
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
    if (debug_adjust) {
      PrintMessageDieOnError("[BWAReadAligner]: Getting qualities for mate", DEBUG);
    }
    // update qualities. For read pairs qual is sum of the two ends' mapq
    int edit;
    int mate_mapq = GetMapq(aln_seq, ref_seq,
			    aligned_seq_quals, &edit);
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
    if (align_debug) {
      PrintMessageDieOnError("[AdjustAlignment]: returning false, " \
                             "refseq index out of range.", DEBUG);
    }
    return false;
  }
  const REFSEQ& refseq = _ref_sequences->at(aligned_read->strid);
  size_t start_pos = !aligned_read->reverse ?
    aligned_read->lStart-1 : aligned_read->rStart;
  // get sequences to pairwise align
  string rseq;
  try {
    rseq = refseq.sequence.substr(start_pos - refseq.start-REFEXTEND, reglen+REFEXTEND);
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
  aligned_read->read_start = start_pos - REFEXTEND;
  aligned_read->read_end = start_pos + reglen;

  // Global alignment read vs. region aligned to
  string aligned_seq_sw, ref_seq_sw;
  int sw_score;
  CIGAR_LIST cigar_list;
  nw(aligned_seq, rseq, aligned_seq_sw, ref_seq_sw,
     false, &sw_score, &cigar_list);

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
  aligned_read->mapq = GetMapq(aligned_seq_sw, ref_seq_sw,
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
    if (align_debug) {
      PrintMessageDieOnError("[AdjustAlignment]: calling GetSTRAllele", DEBUG);
    }
    return GetSTRAllele(aligned_read, cigar_list);
  } catch(std::out_of_range & exception) {
    if (align_debug) {
      PrintMessageDieOnError("[AdjustAlignment]: Problem adjusting read " +
			     aligned_read->ID, WARNING);
    }
    return false;
  }
  if (align_debug) {
    PrintMessageDieOnError("[AdjustAlignment]: returning false, reached " \
                           "end of adjust.", DEBUG);
  }
  return false;
}

int BWAReadAligner::GetMapq(const std::string& aligned_sw_string,
                            const std::string& ref_sw_string,
                            const std::string& aligned_quals,
                            int* edit_dist) {
  *edit_dist = 0;
  size_t qual_index = 0;
  int score = 0;
  if (debug_adjust) {
    PrintMessageDieOnError("[GetMapq]: " + aligned_sw_string, DEBUG);
    PrintMessageDieOnError("[GetMapq]: " + ref_sw_string, DEBUG);
    PrintMessageDieOnError("[GetMapq]: " + aligned_quals, DEBUG);
  }
  for (size_t i = 0; i < aligned_sw_string.length(); i++) {
    const char& alnchar = aligned_sw_string.at(i);
    const char& refchar = ref_sw_string.at(i);
    if (alnchar != '-') {
      if (refchar != '-') {
        // mismatch
        if (alnchar != refchar) {
          (*edit_dist)++;
          score += (static_cast<int>(aligned_quals.at(qual_index))-33);
        }
      } else {
        (*edit_dist)++;
      }
      qual_index++;
    } 
  }
  return score;
}

bool BWAReadAligner::GetSTRAllele(MSReadRecord* aligned_read,
                                  const CIGAR_LIST& cigar_list) {
  if (cigar_debug) {
    PrintMessageDieOnError("[GetSTRAllele] CIGAR " + cigar_list.cigar_string, DEBUG);
  }
  // index where STR starts in the read
  size_t str_index = aligned_read->msStart-aligned_read->read_start + 1;
  // Length of the total STR region
  size_t ms_length = aligned_read->msEnd - aligned_read->msStart;

  // check that not too close to ends
  size_t span = 0;
  for (size_t i = 0; i < cigar_list.cigars.size(); i++) {
    const int& s = cigar_list.cigars.at(i).num;
    const char& t = cigar_list.cigars.at(i).cigar_type;
    if (t == 'M' || t == 'D') span += s;
  }
  size_t str_index_end = aligned_read->read_start + span - aligned_read->msEnd;
  if ((str_index < MIN_DIST_FROM_END || str_index_end < MIN_DIST_FROM_END)) {
    if (debug_adjust) {
      PrintMessageDieOnError("[GetSTRAllele]: STR too close to read end", DEBUG);
    }
    return false;
  }


  // If alignment is too messy, get rid of it
  if (cigar_list.cigars.size() > MAX_CIGAR_SIZE) {
    if (align_debug) {
      PrintMessageDieOnError("[BWAReadAligner]: Discarding read, cigar score too long.", DEBUG);
    }
    return false;
  }

  // same as reference
  if (cigar_list.cigars.size() == 1) {
    if (aligned_read->reverse) {
      aligned_read->detected_ms_nuc =
        reverseComplement(aligned_read->nucleotides).
        substr(str_index - 1, ms_length);
    } else {
      aligned_read->detected_ms_nuc =
        aligned_read->nucleotides.substr(str_index - 1, ms_length);
    }
    aligned_read->diffFromRef = 0;
    return (aligned_read->detected_ms_nuc.length() > MIN_STR_LENGTH);
  }

  // get only cigar score spanning the STR
  const int& str_start_in_cigar =
    aligned_read->msStart - aligned_read->read_start;
  // position into the segment
  int pos = 0;
  // base pairs spanned by this cigar item
  int bp = 0;
  // type of the cigar item
  char cigar_type;
  // index into the cigar score
  size_t cigar_index = 0;
  // diff to go until end of this segment
  int diff = 0;
  // temp cigar list to store when removing flanks
  CIGAR_LIST new_cigar_list;
  // list with only cigars for the STR region
  CIGAR_LIST str_cigar_list;
  // Diff in bp from ref STR
  int diff_from_ref = 0;

  // get rid of left flanking region
  while (pos <= str_start_in_cigar  &&
         cigar_index < cigar_list.cigars.size()) {
    bp = cigar_list.cigars.at(cigar_index).num;
    cigar_type = cigar_list.cigars.at(cigar_index).cigar_type;
    // If match or del, increment position
    if (cigar_type == 'M' || cigar_type == 'D' || cigar_type == 'S') pos += bp;
    // bp to go until we hit STR
    diff = pos - str_start_in_cigar;
    if (diff >= 0) {
      size_t cigar_index_to_include = cigar_index;
      // If left adjacent cigar is not M or S, include it
      if (diff == 0 && (cigar_list.cigars.at(cigar_index).cigar_type == 'M' ||
                        cigar_list.cigars.at(cigar_index).cigar_type == 'S')) {
        cigar_index_to_include += 1;
      } else {
        diff -= cigar_list.cigars.at(cigar_index).num;
      }
      new_cigar_list.cigars.resize(cigar_list.cigars.size() -
                                   cigar_index_to_include);
      copy(cigar_list.cigars.begin() + cigar_index_to_include,
           cigar_list.cigars.end(),
           new_cigar_list.cigars.begin());
      break;
    }
    cigar_index += 1;
  }
  // Update STR cigar taking away left flank
  str_cigar_list.cigars = new_cigar_list.cigars;
  str_cigar_list.ResetString();
  new_cigar_list.cigars.clear();

  // get rid of right flank cigars
  // start at beginning of STR list
  cigar_index = 0;
  // Pos from end of the STR region
  pos = diff;
  int total_str_len = static_cast<int>(ms_length);
  while (pos < total_str_len) {
    if (cigar_index >= str_cigar_list.cigars.size()) {
      return false;
    }
    bp = str_cigar_list.cigars.at(cigar_index).num;
    cigar_type = str_cigar_list.cigars.at(cigar_index).cigar_type;
    if (cigar_type == 'M' || cigar_type == 'D' || cigar_type == 'S')
      pos += bp;
    // Difference between our position and the end of the STR
    diff = pos-total_str_len;
    if (diff >= 0) {
      size_t cigar_index_to_include = cigar_index;
      // If right adjacent is not M or S, include it
      if (cigar_index < str_cigar_list.cigars.size() - 1) {
        const char& next_type = str_cigar_list.cigars.
          at(cigar_index+1).cigar_type;
        if (next_type != 'M' && next_type != 'S' && diff == 0) {
          cigar_index_to_include += 1;
        }
      }
      new_cigar_list.cigars.resize(cigar_index_to_include + 1);
      copy(str_cigar_list.cigars.begin(),
           str_cigar_list.cigars.begin() + cigar_index_to_include + 1,
           new_cigar_list.cigars.begin());
      break;
    }
    cigar_index += 1;
  }
  str_cigar_list.cigars.clear();
  str_cigar_list.cigars = new_cigar_list.cigars;
  str_cigar_list.ResetString();
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
    ms_nuc =  rev_read.substr(str_index - 1, ms_length+diff_from_ref);
  } else {
    ms_nuc =  aligned_read->nucleotides.
      substr(str_index - 1, ms_length+diff_from_ref);
  }
  if (ms_nuc.length() <= MIN_STR_LENGTH) {
    if (align_debug) {
      PrintMessageDieOnError("[BWAReadAligner]: Discarding: detected STR too short.", DEBUG);
    }
    return false;
  }
  aligned_read->diffFromRef = diff_from_ref;
  aligned_read->detected_ms_nuc = ms_nuc;
  if (align_debug) {
    PrintMessageDieOnError("[GetSTRAllele]: returning true", DEBUG);
  }
  return true;
}


BWAReadAligner::~BWAReadAligner() {}

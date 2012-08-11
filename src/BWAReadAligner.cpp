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
// For partial alignment, clear flanking region
// alignments if maps all over, likely
// repetitive
const size_t MAX_ALLOWED_FLANK = 10;
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
  _default_opts->fnr = -1;  // MATE_FNR;

  cigar_debug = false;
  partial_debug = false;
  stitch_debug = false;
}

bool BWAReadAligner::ProcessReadPair(ReadPair* read_pair) {
  if (debug) {
    cerr << "\n[BWAReadAligner]: processing "
         << read_pair->reads.at(0).ID
         << " motif 1 " << read_pair->reads.at(0).repseq;
    if (read_pair->reads.at(0).paired) {
      cerr << " motif 2 " << read_pair->reads.at(1).repseq << endl;
    }
  }
  
  // Initialize status variables
  read_pair->read1_passed_alignment = false;
  read_pair->read2_passed_alignment = false;
  read_pair->found_unique_alignment = false;
  read_pair->aligned_read_num = -1;

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
    if (align_debug) {
      cerr << "[BWAReadAligner]: read1 passed "
           << read_pair->read1_passed_alignment << endl;
      cerr << "[BWAReadAligner]: read2 passed "
           << read_pair->read2_passed_alignment << endl;
    }
    // Step 2: Determine if unique valid alignment
    if ((read_pair->read1_passed_alignment ||
         read_pair->read2_passed_alignment)) {
      if (align_debug) {
        cerr << "[BWAReadAligner]: checking paired alignment"
             << endl;
      }
      ALIGNMENT matealign;
      size_t num_alignments = 0;
      size_t index_of_hit;

      // Get info for the read that aligned
      if (read_pair->read1_passed_alignment) {
        // if both aligned, and 1 is partial,
        // use non-partial one
        if (read_pair->reads.at(0).partial &&
            read_pair->read2_passed_alignment) {
          read_pair->aligned_read_num = 1;
          num_alignments = good_left_alignments_read2.size();
        } else {
          read_pair->aligned_read_num = 0;
          num_alignments = good_left_alignments_read1.size();
        }
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
          cerr << "[BWAReadAligner]: checking for compatible alignment"
               << " read 1 size " << good_left_alignments_read1.size()
               << " read 2 size " << good_left_alignments_read2.size()
               << endl;
        }
        // find compatible alignment
        size_t index_of_mate;
        if (FindCompatibleAlignment(good_left_alignments_read1,
                                    good_left_alignments_read2,
                                    good_right_alignments_read1,
                                    good_right_alignments_read2,
                                    &index_of_hit, &index_of_mate)) {
          if (align_debug) {
            cerr << "[BWAReadAligner]: found compatible alignments " << endl;
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
        if (align_debug) {
          cerr << "[BWAReadAligner]: checking for mate alignment"
               << endl;
        }
        vector<ALIGNMENT> mate_alignments;
        // Trim more harshly here
        string trim_nucs;
        string trim_quals;
        TrimRead(read_pair->reads.at(1-read_pair->aligned_read_num).orig_nucleotides,
                 read_pair->reads.at(1-read_pair->aligned_read_num).orig_qual,
                 &trim_nucs, &trim_quals, MATE_TRIM_QUAL);
        read_pair->reads.at(1-read_pair->aligned_read_num).orig_nucleotides = trim_nucs;
        read_pair->reads.at(1-read_pair->aligned_read_num).orig_qual = trim_quals;
        if (!AlignMate(*read_pair, &mate_alignments,
                       read_pair->reads.at(read_pair->aligned_read_num).repseq)) {
          if (align_debug) {
            cerr << "[BWAReadAligner]: Mate not aligned." << endl;
          }
        } else {
          for (size_t i = 0; i < num_alignments; i++) {
            const ALIGNMENT& lalign = good_left.at(i);
            const ALIGNMENT& ralign = good_right.at(i);
            if (CheckMateAlignment(mate_alignments, lalign, ralign,
                                   &matealign)) {
              if (!read_pair->found_unique_alignment) {
                index_of_hit = i;
                read_pair->found_unique_alignment = true;
              } else {
                // multiple mapper, more than one good hit
                if (align_debug) {
                  cerr << "[BWAReadAligner]: Discarding: multiple mapper."
                       << endl;
                }
                return false;
              }
            }
          }
        }
      }

      // Step 3: Adjust alignment and output
      if (align_debug) {
        cerr << "[BWAReadAligner]: adjust alignment and output" << endl;
      }
      if (read_pair->found_unique_alignment) {
        ALIGNMENT final_left_alignment =
          good_left.at(index_of_hit);
        ALIGNMENT final_right_alignment =
          good_right.at(index_of_hit);
        // try stitching first
        if (align_debug) {
          cerr << "[BWAReadAligner]: Try stitching" << endl;
        }
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
  } else {  // single end reads
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
        if (align_debug) {
          cerr << "\n[BWAReadAligner]: checking single end alignment." << endl;
        }
        // Output the alignment
        read_pair->read1_passed_alignment = true;
        read_pair->found_unique_alignment = true;
        read_pair->aligned_read_num = 0;
        ALIGNMENT dummy_matealign;
        return OutputAlignment(read_pair, good_left_alignments_read1.front(),
                               good_right_alignments_read1.front(),
                               dummy_matealign, false);
      }
    }
  }
  return false;
}

bool BWAReadAligner::ProcessRead(MSReadRecord* read,
                                 vector<ALIGNMENT>* good_left_alignments,
                                 vector<ALIGNMENT>* good_right_alignments) {
  if (align_debug) {
    cerr << "[ProcessRead]: " << read->ID << " " << read->nucleotides << endl;
  }
  // do the flanks consist of all repeats?
  bool left_all_repeats = false;
  bool right_all_repeats = false;
  read->partial = false;
  read->was_partial = false;

  // Detected motif
  const string& repseq = read->repseq;

  // check that we have that repseq
  if (_bwt_references->find(repseq) == _bwt_references->end()) {
    if (align_debug) {
      cerr << "[ProcessRead]: repseq doesn't exist in ref " << repseq << endl;
    }
    return false;
  }

  // Check if flanks are perfect repeats
  if (IsPerfectRepeat(read->left_flank_nuc, repseq) ||
      IsPerfectRepeat(read->left_flank_nuc, reverseComplement(repseq))) {
    read->left_perfect_repeat = true;
    left_all_repeats = true;
  }
  if (align_debug) {
    cerr << "[ProcessRead]: left all repeats " << left_all_repeats << endl;
  }
  if (IsPerfectRepeat(read->right_flank_nuc, repseq) ||
      IsPerfectRepeat(read->right_flank_nuc, reverseComplement(repseq))) {
    read->right_perfect_repeat = true;
    right_all_repeats = true;
  }
  if (align_debug) {
    cerr << "[ProcessRead]: right all repeats " << right_all_repeats << endl;
  }
  // don't continue if both are fully repetitive
  if (left_all_repeats && right_all_repeats) return false;

  if (align_debug) {
    cerr << "[ProcessRead]: align flanks" << endl;
  }
  // Align the flanking regions
  bwa_seq_t* seqs = BWAAlignFlanks(*read);
  bwa_seq_t* seq_left = &seqs[0];
  bwa_seq_t* seq_right = &seqs[1];

  if (align_debug) {
    cerr << "[ProcessRead]: fill in coords" << endl;
  }
  // fill in alignment coordinates
  vector<ALIGNMENT> left_alignments;
  vector<ALIGNMENT> right_alignments;
  if (!left_all_repeats) {
    if (!GetAlignmentCoordinates(seq_left, repseq, &left_alignments)) {
      bwa_free_read_seq(2, seqs);
      if (align_debug) {
        cerr << "[ProcessRead]: no left alignment" << endl;
      }
      return false;
    }
  }
  if (!right_all_repeats) {
    if (!GetAlignmentCoordinates(seq_right, repseq, &right_alignments)) {
      bwa_free_read_seq(2, seqs);
      if (align_debug) {
        cerr << "[ProcessRead]: no right alignment" << endl;
      }
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
    if (align_debug) {
      cerr << "[ProcessRead]: checking partial "
           << left_alignments.size() << " "
           << right_alignments.size() << endl;
    }
    // check if one flanking region aligns uniquely
    if (left_alignments.size() + right_alignments.size() == 1) {
      read->partial = true;
      read->was_partial = true;
      // Case 1: anchored the left side
      if (left_alignments.size() != 0) {
        if (align_debug) {
          cerr << "[ProcessRead]: anchored left" << endl;
        }
        right_all_repeats = true;
        // left side is anchored
        const ALIGNMENT& left_refid =
          left_alignments.front();
        // update right side
        ALIGNMENT right_refid;
        UpdatePartialFlank(left_refid, &right_refid,
                           static_cast<int>(read->nucleotides.length()),
                           static_cast<int>(read->left_flank_nuc.length()),
                           static_cast<int>(read->right_flank_nuc.length()));
        left_refids.push_back(left_refid);
        right_refids.push_back(right_refid);
      } else if (right_alignments.size() != 0) {
        // Case 2: anchored the right side
        if (align_debug) {
          cerr << "[ProcessRead]: anchored right" << endl;
        }
        left_all_repeats = true;
        // Right side is anchored
        const ALIGNMENT& right_refid =
          right_alignments.front();
        // update fields in left
        ALIGNMENT left_refid;
        UpdatePartialFlank(right_refid, &left_refid,
                           static_cast<int>(read->nucleotides.length()),
                           static_cast<int>(read->right_flank_nuc.length()),
                           static_cast<int>(read->left_flank_nuc.length()));
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

  if (align_debug) {
    cerr << "[ProcessRead]: aln sizes "
         << good_left_alignments->size() << " "
         << good_right_alignments->size() << endl;
  }
  // Return true if at least one good alignment
  return (good_left_alignments->size() >= 1 &&
          good_right_alignments->size() >= 1);
}

bwa_seq_t* BWAReadAligner::BWAAlignFlanks(const MSReadRecord& read) {
  if (align_debug) {
    cerr << "[BWAAlignFlanks]: set up" << endl;
  }
  // set up bwa
  int is_comp = _opts->mode&BWA_MODE_COMPREAD;
  bwa_seq_t *seqs, *seq_left, *seq_right;
  seqs = reinterpret_cast<bwa_seq_t*>(calloc(2, sizeof(bwa_seq_t)));
  seq_left = &seqs[0];
  seq_right= &seqs[1];

  if (align_debug) {
    cerr << "[BWAAlignFlanks]: reverse" << endl;
  }
  // left in for historical reasons since BWA reverses sequence
  const string& left_flank_nuc = reverseComplement(read.left_flank_nuc);
  const string& right_flank_nuc = reverseComplement(read.right_flank_nuc);

  if (align_debug) {
    cerr << "[BWAAlignFlanks]: left flank" << endl;
  }
  // LEFT FLANKING REGION
  const string& left_qual =  reverse(read.quality_scores.
                                     substr(0, left_flank_nuc.length()));
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
    cerr << "[BWAAlignFlanks]: right flank" << endl;
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

void BWAReadAligner::UpdatePartialFlank(const ALIGNMENT& aligned_flank,
                                        ALIGNMENT* unaligned_flank,
                                        int nuc_size,
                                        int aligned_flank_size,
                                        int unaligned_flank_size) {
  unaligned_flank->left = !aligned_flank.left;
  unaligned_flank->chrom = aligned_flank.chrom;
  unaligned_flank->start = aligned_flank.start - extend;
  unaligned_flank->end = aligned_flank.end - extend;
  unaligned_flank->id = aligned_flank.id;
  unaligned_flank->copynum = aligned_flank.copynum;
  unaligned_flank->repeat = aligned_flank.repeat;
  unaligned_flank->strand = aligned_flank.strand;
  if (!unaligned_flank->left) {
    unaligned_flank->pos = aligned_flank.pos +
      nuc_size - unaligned_flank_size;
  } else {
    unaligned_flank->pos = aligned_flank.pos +
      + aligned_flank_size - nuc_size;
  }
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
    if (align_debug) {
      cerr << "[AlignMate]: check mate found 0 alignments "
           << nucs << " " << repseq << "\n"
           << qual << endl;
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
      cerr << "[CheckMateAlignment]: check mate finds "
           << it->chrom << " "
           << (abs(it->pos-left_alignment.pos)) << " "
           << it->strand << endl;
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

bool BWAReadAligner::FindCompatibleAlignment(const vector<ALIGNMENT>&
                                             good_left1,
                                             const vector<ALIGNMENT>&
                                             good_left2,
                                             const vector<ALIGNMENT>&
                                             good_right1,
                                             const vector<ALIGNMENT>&
                                             good_right2,
                                             size_t* index_of_hit,
                                             size_t* index_of_mate) {
  if (align_debug) {
    cerr << "[FindCompatibleAlignment]: " \
      "Looking for compatible alignment" << endl;
  }
  bool found_unique = false;
  // For each good left 1, check good left 2
  for (size_t i1 = 0; i1 < good_left1.size(); i1++) {
    for (size_t i2 = 0; i2 < good_left2.size(); i2++) {
      const ALIGNMENT& l1 = good_left1.at(i1);
      const ALIGNMENT& l2 = good_left2.at(i2);
      if (align_debug) {
        cerr << "[FindCompatibleAlignment]: checking pair" << endl;
        cerr << "[FindCompatibleAlignment]: " << l1.chrom << endl;
        cerr << "[FindCompatibleAlignment]: " << l2.chrom << endl;
        cerr << "[FindCompatibleAlignment]: " << l1.pos << endl;
        cerr << "[FindCompatibleAlignment]: " << l2.pos << endl;
        cerr << "[FindCompatibleAlignment]: " << l1.strand << endl;
        cerr << "[FindCompatibleAlignment]: " << l2.strand << endl;
        cerr << "[FindCompatibleAlignment]: " << l1.left << endl;
        cerr << "[FindCompatibleAlignment]: " << l2.left << endl;
      }
      if ((abs(l1.pos-l2.pos) <= MAX_PAIRED_DIFF) &&
          (l1.strand != l2.strand) &&
          (l1.chrom == l2.chrom)) {
        if (found_unique) {
          if (align_debug) {
            cerr << "[BWAReadAligner]: Multiple mapper" << endl;
          }
          return false;
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
      cerr << "[StitchReads]: stitching reads" << endl;
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
      cerr << "[StitchReads]: seq1 " << seq1 << " seq2 " << seq2 << endl;
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
        if (stitch_debug) {
          cerr << "[StitchReads]: adding score " << max_score
               << " " << max_score_index << endl;
        }
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
        if (stitch_debug) {
          cerr << "[StitchReads]: adding score " << max_score
               << " " << max_score_index << endl;
        }
      }
      scores.push_back(score/overlap_len);
    }
    if (stitch_debug) {
      cerr << "[StitchReads]: checking for too many matches" << endl;
    }

    // Check if too many matches
    for (size_t i = 0; i < scores.size(); i++) {
      if ((max_score - scores.at(i) <= STITCH_DIFF) && i != max_score_index+1) {
        if (stitch_debug) {
          cerr << "[StitchReads]: Returning false, too many matches i "
               << i << endl;
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
      if (stitch_debug) {
        cerr << "[StitchReads]: best stitch is backwards. New index "
             << max_score_index << endl;
      }
    }

    if (stitch_debug) {
      cerr << "[StitchReads]: Checking that stitch is good enough " << endl;
    }

    // Check if stitch is good enough
    overlap_len = seq1.length() - max_score_index - 1;
    if ((overlap_len < MIN_STITCH_OVERLAP) ||
        (max_score < STITCH_REQUIRED_SCORE)) {
      if (stitch_debug) {
        cerr << "[StitchReads]: Returning false, score too low "
             << overlap_len << " " << max_score << endl;
      }
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
      if (stitch_debug) {
        cerr << "[StitchReads]: " << na << " " << nb
             << " " << qa << " " << qb << endl;
      }
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
      cerr << "[StitchReads]: orig string 1 " << seq1 << endl;
      cerr << "[StitchReads]: orig string 2 " << seq2 << endl;
      cerr << "[StitchReads]: orig qual 1 " << seq1_qual << endl;
      cerr << "[StitchReads]: orig qual 2 " << seq2_qual << endl;
      cerr << "[StitchReads]: stitched str  " << stitched_string << endl;
    }
    // put stitched info in aligned read
    if (stitch_debug || align_debug) {
      cerr << "[StitchReads]: " << stitched_string << endl;
      cerr << "[StitchReads]: " << stitched_qual << endl;
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
    // If stitch across STR, unset partial
    if (best_stitch_is_backwards) {
      // seq1 is unaligned read
      if (read_pair->reads.at(num_aligned_read).left_all_repeats &&
          !read_pair->reads.at(1-num_aligned_read).left_all_repeats) {
        read_pair->reads.at(num_aligned_read).partial = false;
      }
    } else {
      // seq1 is aligned read
      if (read_pair->reads.at(num_aligned_read).right_all_repeats &&
          !read_pair->reads.at(1-num_aligned_read).right_all_repeats) {
        read_pair->reads.at(num_aligned_read).partial = false;
      }
    }
    return true;
  } catch(std::out_of_range & exception) {
    if (align_debug) {
      cerr << "[StitchReads]: stitching failed. Substring out of range error."
           << endl;
    }
    return false;
  }
}

bool BWAReadAligner::OutputAlignment(ReadPair* read_pair,
                                     const ALIGNMENT& left_alignment,
                                     const ALIGNMENT& right_alignment,
                                     const ALIGNMENT& mate_alignment,
                                     bool treat_as_paired) {
  if (align_debug) {
    cerr << "[BWAReadAligner]: Output alignment" << endl;
  }
  const int& aligned_read_num = read_pair->aligned_read_num;
  read_pair->treat_as_paired = treat_as_paired;

  if (align_debug) {
    cerr << "[BWAReadAligner]: Set info for aligned read ("
         << aligned_read_num << ")" << endl;
  }
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

  if (align_debug) {
    cerr << "[BWAReadAligner]: Checkalignment" << endl;
  }

  // checks to make sure the alignment is reasonable
  // coords make sense on forward strand
  if (!read_pair->reads.at(aligned_read_num).reverse &
      !read_pair->reads.at(aligned_read_num).partial
      & ((read_pair->reads.at(aligned_read_num).lStart >=
          read_pair->reads.at(aligned_read_num).msStart)
         | (read_pair->reads.at(aligned_read_num).rEnd <=
            read_pair->reads.at(aligned_read_num).msEnd))) {
    if (align_debug) {
      cerr << "[BWAReadAligner]: Discarding: forward " \
        "strand coords don't make sense." << endl;
    }
    return false;
  }
  // coords make sense on reverse strand
  if (read_pair->reads.at(aligned_read_num).reverse &
      !read_pair->reads.at(aligned_read_num).partial
      & ((read_pair->reads.at(aligned_read_num).rStart >=
          read_pair->reads.at(aligned_read_num).msStart)
         | (read_pair->reads.at(aligned_read_num).lEnd <=
      read_pair->reads.at(aligned_read_num).msEnd))) {
    if (align_debug) {
      cerr << "[BWAReadAligner]: Discarding: reverse " \
        "strand coords don't make sense." << endl;
    }
    return false;
  }

  // coords of flanks are positive
  if ((read_pair->reads.at(aligned_read_num).rStart <0) |
      (read_pair->reads.at(aligned_read_num).lStart <0) |
      (read_pair->reads.at(aligned_read_num).rEnd <0) |
      (read_pair->reads.at(aligned_read_num).lEnd < 0) |
      (read_pair->reads.at(aligned_read_num).msStart < 0)) {
    if (align_debug) {
      cerr << "[BWAReadAligner]: Discarding: negative coords found."
           << endl;
    }
    return false;
  }

  if (align_debug) {
    cerr << "[BWAReadAligner]: Adjusting alignment" << endl;
  }
  
  // Adjust alignment and STR call
  try {
    if (!AdjustAlignment(&read_pair->reads.at(aligned_read_num),
                         read_pair->reads.
                         at(aligned_read_num).partial,
                         !read_pair->reads.
                         at(aligned_read_num).left_all_repeats,
                         !read_pair->reads.
                         at(aligned_read_num).right_all_repeats)) {
      if (align_debug) {
        cerr << "[BWAReadAligner]: Returning false: "\
          " AdjustAlignment failed." << endl;
      }
      return false;
    }
  } catch(std::out_of_range & exception) {
    if (align_debug) {
      cerr << "[BWAReadAligner]: Returning false: " \
        " AdjustAlignment exception." << endl;
    }
    return false;
  }
  if (align_debug) {
    cerr << "[BWAReadAligner]: Checking unit requirements" << endl;
  }
  // Make sure alignment meets requirements
  if (unit && !read_pair->reads.at(aligned_read_num).partial) {
    if (read_pair->reads.at(aligned_read_num).diffFromRef %
        read_pair->reads.at(aligned_read_num).ms_repeat_best_period != 0)
      if (align_debug) {
        cerr << "[BWAReadAligner]: returning false (unit failed)" << endl;
      }
      return false;
  }
  if (align_debug) {
    cerr << "[BWAReadAligner]: Checking diff and mapq requirements" << endl;
  }
  if ((((abs(read_pair->reads.at(aligned_read_num).diffFromRef) >
        max_diff_ref) ||
       (read_pair->reads.at(aligned_read_num).mapq >= max_mapq)) &&
       !read_pair->reads.at(aligned_read_num).partial)) {
    if (align_debug) {
      cerr << "[BWAReadAligner]: returning false "
           << "(maxdiffref or mapq fail)" << endl;
    }
    return false;
  }
  if (align_debug) {
    cerr << "[BWAReadAligner]: Checking treat as paired" << endl;
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
      cerr << "[BWAReadAligner]: Checking get CIGAR" << endl;
    }
    // get cigar
    try {
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
      const string& rseq = refseq.sequence.
        substr(start_pos - refseq.start, reglen);
      const string& aseq = read_pair->reads.at(1-aligned_read_num).reverse ?
        reverseComplement(read_pair->reads.at(1-aligned_read_num).nucleotides) :
        read_pair->reads.at(1-aligned_read_num).nucleotides;
      const string& aligned_seq_quals =
        read_pair->reads.at(1-aligned_read_num).reverse ?
        reverse(read_pair->reads.at(1-aligned_read_num).quality_scores) :
        read_pair->reads.at(1-aligned_read_num).quality_scores;
      nw(aseq, rseq, aln_seq, ref_seq, false, &score, &cigar_list);
      if (debug_adjust) {
        cerr << "[BWAReadAligner]: Getting qualities for mate" << endl;
      }
      // update qualities. For read pairs qual is sum of the two ends' mapq
      int mate_mapq = GetMapq(aln_seq, ref_seq,
                              aligned_seq_quals);
      const int& read_mapq = read_pair->reads.at(aligned_read_num).mapq;
      read_pair->reads.at(1-aligned_read_num).mapq = mate_mapq+read_mapq;
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
    } catch(std::out_of_range & exception) {
      if (align_debug || debug_adjust) {
        cerr << "[BWAReadAligner]: returning false, " \
          "problem aligning mate." << endl;
      }
      return false;
    }
  }
  return true;
}

bool BWAReadAligner::AdjustAlignment(MSReadRecord* aligned_read,
                                     bool partial,
                                     bool left_aligned,
                                     bool right_aligned) {
  // get reference sequence
  const size_t& reglen = !aligned_read->reverse ?
    (aligned_read->rEnd - aligned_read->lStart) :
    (aligned_read->lEnd - aligned_read->rStart);
  if ((_ref_sequences->find(aligned_read->strid) ==
       _ref_sequences->end())) {
    if (align_debug) {
      cerr << "[AdjustAlignment]: returning false, " \
        "refseq index out of range." << endl;
    }
    return false;
  }
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
  const string& rseq = refseq.sequence.substr(start_pos - refseq.start-REFEXTEND, reglen+REFEXTEND);
  const string& aligned_seq = !aligned_read->reverse ?
    aligned_read->nucleotides :
    reverseComplement(aligned_read->nucleotides);
  const string& aligned_seq_quals = !aligned_read->reverse ?
    aligned_read->quality_scores :
    reverse(aligned_read->quality_scores);

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
  aligned_read->mapq = GetMapq(aligned_seq_sw, ref_seq_sw,
                               aligned_seq_quals);
  // make sure CIGAR is valid
  bool added_s;
  bool cigar_had_s;
  GenerateCorrectCigar(&cigar_list, aligned_read->nucleotides, &added_s, &cigar_had_s);
  aligned_read->cigar = cigar_list.cigars;
  aligned_read->cigar_string = cigar_list.cigar_string;

  if (debug_adjust) {
    cerr << "[AdjustAlignment]: " << aligned_read->ID << endl;
    cerr << "[AdjustAlignment]: " << aligned_seq_sw << endl;
    cerr << "[AdjustAlignment]: " << ref_seq_sw << endl;
    cerr << "[AdjustAlignment]: " << sw_score << endl;
    cerr << "[AdjustAlignment]: " << cigar_list.cigar_string << endl;
  }

  if (align_debug) {
    cerr << "[BWAReadAligner]: checking nwscore" << endl;
  }
  // Check if alignment is reasonably good
  if (sw_score < min_sw_score) {
    if (align_debug) {
      cerr << "[BWAReadAligner]: Discarding read, " \
        "local realignment failed, low sw score." << endl;
    }
    return false;
  }

  if (align_debug) {
    cerr << "[BWAReadAligner]: adjusting partial" << endl;
  }
  // Readjust if partial
  if (partial) {
    // adjust and determine if actually partial
    try {
      if (!AdjustPartialAlignment(aligned_read, cigar_list,
                                  left_aligned, right_aligned,
                                  start_pos, reglen)) {
        if (align_debug) {
          cerr << "[BWAReadAligner]: returning false, " \
            "AdjustPartialAlignment failed." << endl;
        }
        return false;
      }
    } catch(std::out_of_range & exception) {
      if (align_debug) {
        cerr << "[BWAReadAligner]: returning false, " \
          "AdjustPartialAlignment exception." << endl;
      }
      return false;
    }
  }
  // Update STR allele if not partial
  try {
    if (!aligned_read->partial) {
      if (align_debug) {
        cerr << "[AdjustAlignment]: calling GetSTRAllele" << endl;
      }
      return GetSTRAllele(aligned_read, cigar_list);
    } else {
      return true;
    }
  } catch(std::out_of_range & exception) {
    if (align_debug) {
      cerr << "[AdjustAlignment]: returning false, out " \
        "range exception." << endl;
    }
    return false;
  }
  if (align_debug) {
    cerr << "[AdjustAlignment]: returning false, reached "  \
      "end of adjust." << endl;
  }
  return false;
}

int BWAReadAligner::GetMapq(const std::string& aligned_sw_string,
                            const std::string& ref_sw_string,
                            const std::string& aligned_quals) {
  size_t qual_index = 0;
  int score = 0;
  if (debug_adjust) {
    cerr << "[GetMapq]: " << aligned_sw_string << endl;
    cerr << "[GetMapq]: " << ref_sw_string << endl;
    cerr << "[GetMapq]: " << aligned_quals << endl;
  }
  for (size_t i = 0; i < aligned_sw_string.length(); i++) {
    const char& alnchar = aligned_sw_string.at(i);
    const char& refchar = ref_sw_string.at(i);
    if (alnchar != '-') {
      if (refchar != '-') {
        // mismatch
        if (alnchar != refchar) {
          if (debug_adjust) {
            cerr << "[GetMapq]: " << alnchar << " " << refchar << " "
                 << (aligned_quals.at(qual_index)) << endl;
          }
          score += (static_cast<int>(aligned_quals.at(qual_index))-33);
        }
      }
      qual_index++;
    }
  }
  if (debug_adjust) {
    cerr << "[GetMapq]: " << score << endl;
  }
  return score;
}

/*
  Note, this function is quite sketchy due to issues running NW
  on partially spanning alignments. Outputs of this should be treated
  as *rough estimates* that will be corrected in the allelotyping step.
 */
bool BWAReadAligner::AdjustPartialAlignment(MSReadRecord* aligned_read,
                                            const CIGAR_LIST& cigar_list,
                                            bool left_aligned,
                                            bool right_aligned,
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
    if ((cigar_list.cigars.at(i).cigar_type == 'M') ||
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

  if (partial_debug) {
    cerr << "read start " << start_pos << " diff from ref "
         << diff_from_ref << " span " << span << endl;
  }

  // Check if actually partially spanning
  // Case 1: starts at left. see if total span is > (str index+ms len)
  if ((!aligned_read->reverse & left_aligned) |
      (aligned_read->reverse & right_aligned)) {
    strbp = span - (aligned_read->msStart - aligned_read->read_start);
    if (partial_debug) {
      cerr << "case 1, left aligned strbp " << strbp
           << " diff from ref " << diff_from_ref << endl;
    }
    // Check if end extends beyond the read
    if (aligned_read->msEnd - aligned_read->read_start >= span) {
      // Check if reasonable
      if (strbp+diff_from_ref < 0 ||
          strbp+diff_from_ref >=
          static_cast<int>(aligned_read->nucleotides.length()) ||
          strbp < static_cast<int>(MIN_STR_LENGTH)) {
        if (align_debug) {
          cerr << "[BWAReadAligner]: Discarding read, partial " \
            "alignment doesn't make sense." << endl;
        }
        return false;
      }
      if (partial_debug) {
        cerr << "In case 1, updating read data" << endl;
        cerr << "strbp " << strbp << " diff " << diff_from_ref << endl;
      }
      // Update read data
      aligned_read->partial = true;
      aligned_read->was_partial = true;
      aligned_read->diffFromRef = (strbp+diff_from_ref) - ms_length;
      aligned_read->detected_ms_nuc = !aligned_read->reverse ?
        aligned_read->nucleotides.
        substr(aligned_read->nucleotides.length() -
               strbp-diff_from_ref, strbp+diff_from_ref) :
        reverseComplement(aligned_read->nucleotides).
        substr(aligned_read->nucleotides.length() -
               strbp-diff_from_ref, strbp+diff_from_ref);
      return true;
    } else if (!aligned_read->left_perfect_repeat &&
               !aligned_read->right_perfect_repeat) {
      if (debug_adjust) {
        cerr << "[AdjustPartialAlignment]: actually not partial"
             << endl;
      }
      // aligned_read->partial = false;
      aligned_read->diffFromRef = (strbp+diff_from_ref) - ms_length;
      return true;
    } else if (aligned_read->right_perfect_repeat ||
               aligned_read->left_perfect_repeat) {
      // Update diff from ref, know alignment was bad
      aligned_read->diffFromRef = (aligned_read->nucleotides.length() -
                                   (aligned_read->msStart - aligned_read->read_start)
                                   - ms_length);
      return true;
    }
  }

  // Case 2: anchored at the right
  if ((!aligned_read->reverse & right_aligned) |
      (aligned_read->reverse & left_aligned)) {
    strbp = span - (start_pos+reglen-aligned_read->msEnd);
    if (partial_debug) {
      cerr << "case 2, right aligned strbp " << strbp << " span "
           << span << " reglen " << reglen << endl;
    }
    // Check if beginning extends beyond the msStart
    if ((start_pos  + reglen - aligned_read->msStart) >= span) {
      // check if reasonable
      if (strbp+diff_from_ref < 0 ||
          strbp+diff_from_ref >=
          static_cast<int>(aligned_read->nucleotides.length()) ||
          strbp <= static_cast<int>(MIN_STR_LENGTH)) {
        if (align_debug) {
          cerr << "[BWAReadAligner]: Discarding read, partial " \
            "alignment doesn't make sense." << endl;
        }
        return false;
      }
      if (partial_debug) {
        cerr << "In case 2, updating read data" << endl;
        cerr << "strbp " << strbp << " diff " << diff_from_ref << endl;
      }
      // update read info
      aligned_read->partial = true;
      aligned_read->was_partial = true;
      aligned_read->diffFromRef = (strbp +diff_from_ref)-ms_length;
      aligned_read->detected_ms_nuc = aligned_read->reverse ?
        aligned_read->nucleotides.substr(0, strbp+diff_from_ref) :
        reverseComplement(aligned_read->nucleotides).
        substr(0, strbp+diff_from_ref);
      return true;
    } else if (!aligned_read->left_perfect_repeat &&
               !aligned_read->right_perfect_repeat) {
      if (debug_adjust) {
        cerr << "[AdjustPartialAlignment]: actually not partial" << endl;
      }
      aligned_read->diffFromRef = (strbp+diff_from_ref) - ms_length;
      // aligned_read->partial = false;
      // Need to put something in detected ms nuc. In the partial
      // case this field is unreliable.
      aligned_read->detected_ms_nuc = aligned_read->reverse ?
        aligned_read->nucleotides : reverseComplement(aligned_read->nucleotides);
      return true;
    } else if (aligned_read->left_perfect_repeat ||
               aligned_read->right_perfect_repeat ) {
      // Update diff from ref, know alignment was bad
      // Need to put something in detected ms nuc. In the partial
      // case this field is unreliable.
      aligned_read->detected_ms_nuc = aligned_read->reverse ?
        aligned_read->nucleotides : reverseComplement(aligned_read->nucleotides);
      aligned_read->diffFromRef = (start_pos+reglen-aligned_read->msEnd)-ms_length;
      return true;
    }
  }
  if (debug_adjust) {
    cerr << "[AdjustPartialAlignment]: left perfect ? "
         << aligned_read->left_perfect_repeat << endl;
    cerr << "[AdjustPartialAlignment]: right perfect ? "
         << aligned_read->right_perfect_repeat << endl;
    cerr << "[AdjustPartialAlignment]: end of adjust partial" << endl;
  }
  return false;
}

bool BWAReadAligner::GetSTRAllele(MSReadRecord* aligned_read,
                                  const CIGAR_LIST& cigar_list) {
  if (cigar_debug) {
    cerr << "CIGAR " << cigar_list.cigar_string << endl;
  }
  // index where STR starts in the read
  size_t str_index = aligned_read->msStart-aligned_read->read_start + 1;
  // Length of the total STR region
  size_t ms_length = aligned_read->msEnd - aligned_read->msStart;

  // check that not too close to ends, else call it partial
  size_t span = 0;
  for (size_t i = 0; i < cigar_list.cigars.size(); i++) {
    const int& s = cigar_list.cigars.at(i).num;
    const char& t = cigar_list.cigars.at(i).cigar_type;
    if (t == 'M' || t == 'D') span += s;
  }
  size_t str_index_end = aligned_read->read_start + span - aligned_read->msEnd;
  if (debug_adjust) {
    cerr << "[GetSTRAllele]: dist start " << str_index << " "
         << " dist end " << str_index_end << endl;
  }
  if ((str_index < MIN_DIST_FROM_END || str_index_end < MIN_DIST_FROM_END)
      && aligned_read->was_partial) {
    if (debug_adjust) {
      cerr << "[GetSTRAllele]: changing to partial, STR too close to read end"
           << endl;
    }
    aligned_read->partial = true;
    return true;
  }


  // If alignment is too messy, get rid of it
  if (cigar_list.cigars.size() > MAX_CIGAR_SIZE) {
    if (align_debug) {
      cerr << "[BWAReadAligner]: Discarding read, cigar score too long."
           << endl;
    }
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
    if (align_debug) {
      cerr << "[GetSTRAllele]: str len "
           << (aligned_read->detected_ms_nuc.length())
           << endl;
    }
    return (aligned_read->detected_ms_nuc.length() > MIN_STR_LENGTH);
  }

  // get only cigar score spanning the STR
  const int& str_start_in_cigar =
    aligned_read->msStart - aligned_read->read_start;
  if (cigar_debug) {
    cerr << "ms start " << aligned_read->msStart << endl;
    cerr << "read start " << aligned_read->read_start << endl;
    cerr << "start in cigar " << str_start_in_cigar << endl;
  }
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
    if (cigar_debug) {
      cerr << "before: diff " << diff <<  " pos " << pos
           << " str start " << str_start_in_cigar << endl;
    }
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

  if (cigar_debug) {
    cerr << "after: diff " << diff <<  " pos " << pos
         << " str start " << str_start_in_cigar << endl;
    cerr << "str cigar string minus left flank "
         << str_cigar_list.cigar_string << endl;
  }

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
      if (cigar_debug) {
        cerr << " diff " << diff << " cigar index "
             << cigar_index << " pos " << pos << endl;
      }
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
    if (cigar_debug) {
      cerr << "STR CIGAR after remove right flank "
           << str_cigar_list.cigar_string << endl;
    }

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
  if (cigar_debug) {
    cerr << "diff " << diff_from_ref << endl;
    cerr << "ms len " << ms_length << endl;
    cerr << "str index " << str_index << endl;
    cerr << " nucs " << aligned_read->nucleotides << endl;
  }

  // set STR region
  string ms_nuc;
  if (aligned_read->reverse) {
    string rev_read = reverseComplement(aligned_read->nucleotides);
    ms_nuc =  rev_read.substr(str_index, ms_length+diff_from_ref);
  } else {
    ms_nuc =  aligned_read->nucleotides.
      substr(str_index, ms_length+diff_from_ref);
  }
  if (ms_nuc.length() <= MIN_STR_LENGTH) {
    if (align_debug) {
      cerr << "[BWAReadAligner]: Discarding: detected STR too short." << endl;
    }
    return false;
  }
  aligned_read->diffFromRef = diff_from_ref;
  aligned_read->detected_ms_nuc = ms_nuc;
  if (align_debug) {
    cerr << "[GetSTRAllele]: returning true" << endl;
  }
  return true;
}


BWAReadAligner::~BWAReadAligner() {}

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

#include <stdexcept>

#include "src/AlignmentUtils.h"
#include "src/common.h"
#include "src/runtime_parameters.h"

using namespace std;

// Minimum number of bp for stitch overlap
const size_t MIN_STITCH_OVERLAP = 16;
// Percent identity required to stitch
const float STITCH_REQUIRED_SCORE = 0.8;
// Allowed difference in score between returned stitch
// and next best stitch
const float STITCH_DIFF = 0.1;
// min allowed distance from STR boundary to read ends
size_t MIN_DIST_FROM_END = 8;
// more than this is likely bad alignment
const size_t MAX_CIGAR_SIZE = 5;

namespace AlignmentUtils {
  int GetMapq(const string& aligned_sw_string,
	      const string& ref_sw_string,
	      const string& aligned_quals,
	      int* edit_dist) {
    *edit_dist = 0;
    size_t qual_index = 0;
    int score = 0;
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

  bool StitchReads(ReadPair* read_pair,
		   ALIGNMENT* left_alignment,
		   ALIGNMENT* right_alignment) {
    try { // TODO remove this try/catch
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
      
      // Check if too many matches
      for (size_t i = 0; i < scores.size(); i++) {
	if ((max_score - scores.at(i) <= STITCH_DIFF) && i != max_score_index+1) {
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
      
      // put stitched info in aligned read
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
      PrintMessageDieOnError("Stitching failed " + read_pair->reads.at(0).ID, WARNING);
      return false;
    }
  }

  bool GetSTRAllele(MSReadRecord* aligned_read,
		    const CIGAR_LIST& cigar_list) {
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
      return false;
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
	  substr(str_index - 1, ms_length);
      } else {
	aligned_read->detected_ms_nuc =
	  aligned_read->nucleotides.substr(str_index - 1, ms_length);
      }
      aligned_read->diffFromRef = 0;
      return (aligned_read->detected_ms_nuc.length() >= MIN_STR_LENGTH);
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
      return false;
    }
    aligned_read->diffFromRef = diff_from_ref;
    aligned_read->detected_ms_nuc = ms_nuc;
    return true;
  }
}

/*
Copyright (C) 2014 Thomas Willems <twillems@mit.edu>

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

#include <vector>
#include <string>
#include <sstream>

#include "src/api/BamAlignment.h"
#include "src/AlignmentFilters.h"
#include "src/common.h"
#include "src/ZAlgorithm.h"

using namespace std;

namespace AlignmentFilters {
  template<typename CigarIterator> int GetDistToIndel(CigarIterator iter, CigarIterator end){
    // Process leading clipping ops
    if (iter != end && iter->Type == 'H')
      iter++;
    if (iter != end && iter->Type == 'S')
      iter++;
    
    int dist = 0;
    while (iter != end){
      char type = iter->Type;
      if (type == 'M')
        dist += iter->Length;
      else if (type == 'I' || type == 'D')
        return dist;
      else if (type == 'S' || type == 'H')
        return -1;
      else {
        string msg = "Invalid CIGAR char";
        msg += type;
        PrintMessageDieOnError(msg, ERROR);
      }
      iter++;
    }
    return -1;
  }
  
  pair<int,int> GetEndDistToIndel(AlignedRead* aln){
    vector<BamTools::CigarOp>::iterator cigar_iter;
    vector<BamTools::CigarOp>::iterator cigar_end;
    vector<int> vals;
    int head_dist = GetDistToIndel(aln->cigar_ops.begin(),  aln->cigar_ops.end());
    int tail_dist = GetDistToIndel(aln->cigar_ops.rbegin(), aln->cigar_ops.rend());
    return pair<int,int>(head_dist, tail_dist);
  }
  
  pair<int,int> GetNumEndMatches(AlignedRead* aln, const string& ref_seq, int ref_seq_start){
    if (aln->read_start < ref_seq_start)
      return pair<int,int>(-1,-1);
    
    unsigned int read_index = 0;
    unsigned int ref_index  = aln->read_start-ref_seq_start;
    vector<BamTools::CigarOp>::iterator cigar_iter = aln->cigar_ops.begin();
    bool beginning = true;
    int match_run  = 0;
    int head_match = 0;
    
    // Process leading clip CIGAR types
    if (cigar_iter != aln->cigar_ops.end() && cigar_iter->Type == 'H')
      cigar_iter++;
    if (cigar_iter != aln->cigar_ops.end() && cigar_iter->Type == 'S'){
      read_index += cigar_iter->Length;
      cigar_iter++;
    }
    
    // Process CIGAR items as long as read region lies within reference sequence bounds
    while (cigar_iter != aln->cigar_ops.end() && ref_index < ref_seq.size() && read_index < aln->nucleotides.size()){
      if (cigar_iter->Type == 'M'){
        if (ref_index + cigar_iter->Length > ref_seq.size()) 
          return pair<int,int>(-1, -1);
        if (read_index + cigar_iter->Length > aln->nucleotides.size())
          PrintMessageDieOnError("Nucleotides for aligned read don't correspond to the CIGAR string", ERROR);
        for (unsigned int len = cigar_iter->Length; len > 0; len--){
          if (ref_seq[ref_index] == aln->nucleotides[read_index])
            match_run++;
          else {
            if (beginning) head_match = match_run;
            beginning = false;
            match_run = 0;
          }
          read_index++;
          ref_index++;
        }
      }
      else if (cigar_iter->Type == 'I'){
        if (beginning) head_match = match_run;
        beginning   = false;
        match_run   = 0;
        read_index += cigar_iter->Length;
      }
      else if (cigar_iter->Type == 'D'){
        if (beginning) head_match = match_run;
        beginning  = false;
        match_run  = 0;
        ref_index += cigar_iter->Length;
      }
      else if (cigar_iter->Type == 'S' || cigar_iter->Type == 'H')
        break;
      else {
        string msg = "Invalid CIGAR char";
        msg += cigar_iter->Type;
        PrintMessageDieOnError(msg, ERROR);
      }
      cigar_iter++;
    }
    
    // Process trailing clip CIGAR types
    if (cigar_iter != aln->cigar_ops.end() && cigar_iter->Type == 'S'){
      read_index += cigar_iter->Length;
      cigar_iter++;
    }
    if (cigar_iter != aln->cigar_ops.end() && cigar_iter->Type == 'H')
      cigar_iter++;
    
    // Ensure that we processed all CIGAR options
    if (cigar_iter != aln->cigar_ops.end()){
      if (ref_index >= ref_seq.size())
        return pair<int,int>(-1,-1);
      else
        PrintMessageDieOnError("Improperly formatted CIGAR string", ERROR);
    }
    
    // Ensure that CIGAR string corresponded to aligned bases
    if (read_index != aln->nucleotides.size()){
      if (ref_index >= ref_seq.size())
        return pair<int,int>(-1,-1);
      else
        PrintMessageDieOnError("CIGAR string does not correspond to alignment bases", ERROR);
    }
    
    if (beginning)
      return pair<int,int>(match_run, match_run);
    else
      return pair<int,int>(head_match, match_run);
  } 
  
  
  /* 
     Stores the sequence, start and end position of the read after removing clipped bases
     using the provided references
  */
  void GetUnclippedInfo(AlignedRead* aln, string& bases, int& unclipped_start, int& unclipped_end){
    unclipped_start = aln->read_start;
    unclipped_end   = aln->read_start-1;
    bool begin      = true;
    int start_index = 0, num_bases = 0;
    for(vector<BamTools::CigarOp>::iterator cigar_iter = aln->cigar_ops.begin(); cigar_iter != aln->cigar_ops.end(); cigar_iter++){
      switch(cigar_iter->Type) {
      case 'D':
        unclipped_end += cigar_iter->Length;
        begin          = false;
        break;
      case 'H':
        break;
      case 'S':
        if (begin) start_index += cigar_iter->Length;
        break;
      case 'M':
        unclipped_end += cigar_iter->Length;
        num_bases     += cigar_iter->Length;
        begin          = false;
        break;
      case 'I':
        num_bases += cigar_iter->Length;
        begin      = false;
        break;
      default:
        string msg = "Invalid CIGAR char ";
        msg += cigar_iter->Type;
        PrintMessageDieOnError(msg, ERROR);
        break;
      }
    }
    bases = aln->nucleotides.substr(start_index, num_bases);
  }
    
  bool HasLargestEndMatches(AlignedRead* aln, const string& ref_seq, int ref_seq_start, int max_external, int max_internal){
    // Extract sequence, start and end coordinates of read after clipping
    string bases;
    int start, end;
    GetUnclippedInfo(aln, bases, start, end);
    
    // Check that the prefix match is the longest
    if (start >= ref_seq_start && start < ref_seq_start + static_cast<int>(ref_seq.size())){
      int start_index = start - ref_seq_start;
      int start       = max(0, start_index - max_external);
      int stop        = min(static_cast<int>((ref_seq.size()-1)), start_index + max_internal);
      vector<int> match_counts;
      ZAlgorithm::GetPrefixMatchCounts(bases, ref_seq, start, stop, match_counts);
      
      int align_index = start_index - start;
      int num_matches = match_counts[align_index];
      for (int i = 0; i < static_cast<int>(match_counts.size()); i++){
        if (i == align_index)
          continue;
        if (match_counts[i] >= num_matches)
          return false;
      }
    }
    
    // Check that the suffix match is the longest
    if (end >= ref_seq_start && end < ref_seq_start + static_cast<int>(ref_seq.size())){
      int end_index = end - ref_seq_start;
      int start     = max(0, end_index - max_internal);
      int stop      = min(static_cast<int>(ref_seq.size()-1), end_index + max_external);
      vector<int> match_counts;
      ZAlgorithm::GetSuffixMatchCounts(bases, ref_seq, start, stop, match_counts);
      
      int align_index = end_index - start;
      int num_matches = match_counts[align_index];
      for (int i = 0; i < static_cast<int>(match_counts.size()); i++){
        if (i == align_index)
          continue;
        if (match_counts[i] >= num_matches)
          return false;
      }
    }       
    return true;
  }

  int GetMaxRepeatsInEnds(AlignedRead* aln, size_t bp_from_end) {
    // Get ends
    if (aln->nucleotides.length() < bp_from_end) {
      return 0;
    }
    const std::string left_end = aln->nucleotides.substr(0, bp_from_end);
    const std::string right_end = aln->nucleotides.substr(aln->nucleotides.length() - bp_from_end, bp_from_end);
    // Try each shift of motif
    int left_occurrences = CountOccurrences(left_end, aln->repseq);
    int right_occurrences = CountOccurrences(right_end, aln->repseq);
    for (size_t shift=1; shift<aln->repseq.size(); shift++) {
      const string shifted_motif = aln->repseq.substr(shift, aln->repseq.size() - shift) +
        aln->repseq.substr(0, shift);
      int lo = CountOccurrences(left_end, shifted_motif);
      int ro = CountOccurrences(right_end, shifted_motif);
      if (lo > left_occurrences) {
        left_occurrences = lo;
      }
      if (ro > right_occurrences) {
        right_occurrences = ro;
      }
    }
    if (left_occurrences > right_occurrences) {
      return left_occurrences;
    }
    return right_occurrences;
  }

  void GetDistDiffFromEnd(AlignedRead* aln) {
    // Get span of read
    int span = 0;
    for (vector<BamTools::CigarOp>::const_iterator it = aln->cigar_ops.begin();
         it != aln->cigar_ops.end(); it++) {
      if (it->Type == 'D' || it->Type == 'M' ||
          it->Type == '=' || it->Type == 'X') {
        span += it->Length;
      }
    }
    int left_dist = aln->msStart - aln->read_start;
    int right_dist = aln->read_start + span - aln->msEnd;
    aln->dist_from_end = left_dist - right_dist;
  }
}


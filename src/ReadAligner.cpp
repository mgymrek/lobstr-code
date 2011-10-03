/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <assert.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string>

#include "common.h"
#include "ReadAligner.h"
#include "runtime_parameters.h"

using namespace std;

ReadAligner::ReadAligner(map<string,GST*> *_gstsL, map<string,GST*> *_gstsR,
			 map<int,MSRecord> * _msDict) {
  gstsL = _gstsL;
  gstsR = _gstsR;
  msDict = _msDict;
}

ReadAligner::~ReadAligner() {}

bool ReadAligner::ProcessRead(MSReadRecord* read){
  if (align_debug) {
    cout << "\n processing " << read->nucleotides << endl;
  }
  int aligned_forward;
  int aligned_reverse;
  bool reverse = false;
  
  NodeAlignment left_node_alignment;
  NodeAlignment right_node_alignment;
  NodeAlignment left_node_alignment_reverse;
  NodeAlignment right_node_alignment_reverse;
  int alignment_id;
  int alignment_id_reverse;
  int left_mismatch;
  int right_mismatch;
  int left_mismatch_reverse;
  int right_mismatch_reverse;

  aligned_forward = AlignRead(*read, false, &left_node_alignment,
			      &right_node_alignment, &alignment_id,
			      &left_mismatch, &right_mismatch);
  aligned_reverse = AlignRead(*read, true, &left_node_alignment_reverse,
			      &right_node_alignment_reverse,
			      &alignment_id_reverse, &left_mismatch_reverse,
			      &right_mismatch_reverse);
  if (align_debug) {
    cout << "forward " << aligned_forward << " reverse " << aligned_reverse << endl;
  }
  int right_flank_start;
  int left_flank_start;
  int right_flank_end;
  int left_flank_end;
  int left_flank_length = read->left_flank_nuc.size();
  int right_flank_length = read->right_flank_nuc.size();
  int left_index_from_start = read->left_flank_index_from_start;
  int right_index_from_end = read->right_flank_index_from_end;
  if (aligned_forward + aligned_reverse == 1) {
    if (aligned_reverse) {
      reverse = true;
      left_node_alignment = left_node_alignment_reverse;
      right_node_alignment = right_node_alignment_reverse;
      alignment_id = alignment_id_reverse;
      left_flank_length = read->right_flank_nuc.size();
      right_flank_length = read->left_flank_nuc.size();
      left_mismatch = right_mismatch_reverse;
      right_mismatch = left_mismatch_reverse;
      left_index_from_start = read->right_flank_index_from_end;
      right_index_from_end = read->left_flank_index_from_start;
    }
    right_flank_start = right_node_alignment.GetStartLocation(alignment_id);
    left_flank_start = left_node_alignment.GetStartLocation(alignment_id);
    right_flank_end = right_flank_start + right_flank_length + 1;
    left_flank_end = left_flank_start + left_flank_length;
  } else {
    // didn't find unique alignment
    return false;
  }

  // get the STR locus information
  MSRecord msInfo;
  if (msDict->find(alignment_id) == msDict->end()) {
    return false;
  }
  msInfo = msDict->at(alignment_id);
  read->chrom = msInfo.chrom;
  read->msStart = msInfo.start;
  read->msEnd = msInfo.end;
  read->msRepeat = msInfo.repeat;
  read->refCopyNum = msInfo.copynum;
  
  /*
  if(reverse){
    ref_length_at_locus =
      (right_node_alignment.GetStartLocation(alignment_id) +
	read->left_flank_nuc.length()) -
       left_node_alignment.GetStartLocation(alignment_id);
      //(msInfo.start + right_flank_start + 1 + read->left_flank_nuc.size()) -
      //(msInfo.start - extend + left_flank_start +1) +1;
    read->rEnd = (msInfo.start+right_flank_start  +
		  read->left_flank_nuc.size() + right_index_from_end);
  }
  else{
    ref_length_at_locus =
    (right_node_alignment.GetStartLocation(alignment_id) +
      read->right_flank_nuc.length()) -
     left_node_alignment.GetStartLocation(alignment_id);
    //(msInfo.start + right_flank_start + 1 + read->right_flank_nuc.size())-
    //  (msInfo.start - extend + left_flank_start +1) +1;
    read->rEnd = (msInfo.start+right_flank_start +
		  read->right_flank_nuc.size() + right_index_from_end);
		  }*/

  // Note (msIndo.start - extend) gives 0 index
  read->lStart = msInfo.start - extend + left_flank_start;
  read->lEnd = msInfo.start - extend + left_flank_end;
  read->rStart = msInfo.start - extend + right_flank_start + 1;
  read->rEnd = msInfo.start + right_flank_end - 1;
  read->reverse = reverse;
  read->lDist = left_mismatch;
  read->rDist = right_mismatch;
  read->name = msInfo.name;

  // get length diff from ref
  int read_length_at_locus = read->left_flank_nuc.length() +
    read->detected_ms_region_nuc.length() +
    read->right_flank_nuc.length() - 1;
  
  int ref_length_at_locus = read->rEnd - read->lStart;
  read->diffFromRef = read_length_at_locus - ref_length_at_locus + 1;

  // update the left start and right end for the case where
  // flanking regions were trimmed
  read->lStart = read->lStart - left_index_from_start;
  read->rEnd = read->rEnd + right_index_from_end;

  // If the alignment is good, return it
  if (align_debug) {
    cout << "read aligned after trimming "
	 << read->nucleotides.substr(read->left_flank_index_from_start,read_length_at_locus) << endl;
    cout << "lstart " << read->lStart 
	 << " lend " << read->lEnd
	 << " rstart " << read->rStart
	 << " rend " << read->rEnd << endl;
    cout << "diff from ref " << abs(read->diffFromRef)
	 << " read len at locus " << read_length_at_locus
	 << " ref len at locus " << ref_length_at_locus << endl;}
  if (abs(read->diffFromRef) <= max_diff_ref &&
     (read->msStart - read->lStart > 0)){return true;}
  // added the >0 part because it messes up sam format and is obviously wrong
  return false;
}

bool ReadAligner::getMSSeq(const string& nucs, int k, string* repeat) {
  if (nucs.size() <= k*2) {return false;}
  unsigned const int base = 10;
  unsigned long long xPowOfBase = 1;
  int i = 0;
  for (i = 1; i <=k; ++i) {
    xPowOfBase *= base;
  }
  unsigned long long firstXLengthSubString = 0;
  for (i = 0; i < k; ++i) {
    firstXLengthSubString *=base;
    firstXLengthSubString += nucs[i];
  }

  unsigned long long nextXLengthSubstring = firstXLengthSubString;
  map<unsigned long long, int > hashTable;
  int max_count = 0;
  int first_pos = -1;
  for (; i <=nucs.size(); ++i) {
    if (hashTable.find(nextXLengthSubstring) != hashTable.end()) {
      ++hashTable[nextXLengthSubstring];
    } else {
      hashTable.insert(make_pair(nextXLengthSubstring, 1));
    }
    int cur_count = hashTable[nextXLengthSubstring];
    if (cur_count > max_count) {
      max_count = cur_count;
      first_pos = i-k;
    }
    if (i != nucs.size()) {
	nextXLengthSubstring *=base;
	nextXLengthSubstring += nucs[i];
	nextXLengthSubstring -= nucs[i-k]*xPowOfBase;
    }
  }
  if (first_pos == -1) {
    return false;
  }
  string kmer = nucs.substr(first_pos, k); 
  /*
  map<string,int> countKMers;
  string subseq;
  string kmer;
  int maxkmer = 0;
  subseq.resize(k);
  // TODO can use rabin karp here?
  for(int i = 0; i < nucs.size()-k; i++){
    subseq = nucs.substr(i,k);
    countKMers[subseq]++;
    if (countKMers.at(subseq) > maxkmer){
      kmer = subseq;
      maxkmer = countKMers.at(subseq);
    }
  }*/
  getCanonicalMS(kmer, repeat);
  return true;
}

int ReadAligner::AlignRead(const MSReadRecord& read, bool reverse,
			   NodeAlignment* left_node_alignment,
			   NodeAlignment* right_node_alignment,
			   int* alignment_id, int* left_mismatch,
			   int* right_mismatch) {
  // number of mismatches in alignment
  int mismatch = 0;

  // take care of reverse complementing
  string left_flank_nucs;
  string right_flank_nucs;
  string repeat_nucs;
  left_flank_nucs = read.left_flank_nuc;
  right_flank_nucs = read.right_flank_nuc;
  repeat_nucs = read.detected_ms_region_nuc;
  if (reverse) {
    // reverse everything
    left_flank_nucs = reverseComplement(read.right_flank_nuc);
    right_flank_nucs = reverseComplement(read.left_flank_nuc);
    repeat_nucs = reverseComplement(read.detected_ms_region_nuc);
  }
  // determine the canonical STR sequence
  string str_repeat;
  if (!getMSSeq(repeat_nucs, read.ms_repeat_best_period,&str_repeat)) { 
    return false;
  }
  if (align_debug) { cout << "msseq " << str_repeat << endl;}
  // align left flank
  if (IsPerfectRepeat(left_flank_nucs, str_repeat) ||
      IsPerfectRepeat(right_flank_nucs, str_repeat)) {
    return false;
  }
  list<NodeAlignment> left_alignments;
  if (gstsL->find(str_repeat) !=
      gstsL->end()) {
    while(mismatch <= allowed_mismatches && left_alignments.size() == 0) {
      if (align_debug) {
	cout << "align left " << left_flank_nucs << " " << str_repeat << endl;
      }
      if (left_flank_nucs.length() < min_length_to_allow_mismatches) {
	left_alignments = gstsL->at(str_repeat)->
	  fuzzyMatch(left_flank_nucs, 0, max_align, 0);
	mismatch++;
	break;
      } else {
	left_alignments = gstsL->at(str_repeat)->
	  fuzzyMatch(left_flank_nucs, mismatch, max_align, 0);
	mismatch++;
      }
    }
  } else { 
    return 0;
  }
  if (left_alignments.size() == 0) {
    return 0;
  }
  *left_mismatch = mismatch - 1;
  mismatch--; // counted one too many on last loop iteration above

  // align right flank
  int right_allowed_mismatches = allowed_mismatches - mismatch;
  mismatch = 0;
  list<NodeAlignment> right_alignments;
  if (gstsR->find(str_repeat) !=
      gstsR->end()) {
    while(mismatch <= right_allowed_mismatches && right_alignments.size() == 0) {
      if (right_flank_nucs.length() < min_length_to_allow_mismatches) {
	right_alignments = gstsR->at(str_repeat)->
	  fuzzyMatch(right_flank_nucs, 0, max_align, 0);
	mismatch++;
	break;
      } else {
	right_alignments = gstsR->at(str_repeat)->
	  fuzzyMatch(right_flank_nucs, mismatch,
		     max_align, 0);
	mismatch++;
      }
    }
  } else { 
    return 0;
  }
  if (right_alignments.size() == 0) {
    return 0;
  }
  *right_mismatch = mismatch - 1;

  // see if we found a unique alignment
  // get unique ids in the left and right flanks

  map<int, NodeAlignment*> left_unique_identifiers;
  GetUniqueIdentifiers(&left_alignments, &left_unique_identifiers);
  if (left_unique_identifiers.size() == 0) {return 0;}
  map<int, NodeAlignment*> right_unique_identifiers;
  GetUniqueIdentifiers(&right_alignments, &right_unique_identifiers);
  if (right_unique_identifiers.size() == 0) {return 0;}

  // get intersection of identifiers in left and right
  list<int> shared_identifiers;
  GetSharedKeys(left_unique_identifiers, right_unique_identifiers,
		&shared_identifiers);

  if (shared_identifiers.size() != 1) {
    return 0;
  }
  *alignment_id = shared_identifiers.front();
  assert(left_unique_identifiers.find(*alignment_id) !=
	 left_unique_identifiers.end());
  assert(right_unique_identifiers.find(*alignment_id) !=
	 right_unique_identifiers.end());
  *left_node_alignment = *left_unique_identifiers.at(*alignment_id);
  *right_node_alignment = *right_unique_identifiers.at(*alignment_id);
  assert(left_node_alignment->string_id_to_alignment.find(*alignment_id) !=
	 left_node_alignment->string_id_to_alignment.end());
  assert(right_node_alignment->string_id_to_alignment.find(*alignment_id) !=
	 right_node_alignment->string_id_to_alignment.end());
  return 1;
}

void ReadAligner::GetUniqueIdentifiers(list<NodeAlignment>* node_alignments,
				       map<int, NodeAlignment*>* unique_identifiers) {
  for (list<NodeAlignment>::iterator it = node_alignments->begin();
       it != node_alignments->end(); ++it) {
    list<int> string_identifiers_list;
    (*it).GetStringIdentifiers(&string_identifiers_list);
    for (list<int>::iterator ids_it = string_identifiers_list.begin();
	 ids_it != string_identifiers_list.end(); ++ids_it) {
      // check if it's in the map, if it is remove what's there
      // and don't add anything, we don't want any duplicates
      if (unique_identifiers->find(*ids_it) != unique_identifiers->end()) {
	unique_identifiers->erase(*ids_it);
	continue;
      } 
      // this entry isn't there yet
      unique_identifiers->insert(pair<int, NodeAlignment*>(*ids_it, &(*it)));
    }
  }
}

void ReadAligner::GetSharedKeys(const map<int, NodeAlignment*>& map1,
				const map<int, NodeAlignment*>& map2,
				list<int>* shared_keys) {
  list<int> keys1;
  for (map<int, NodeAlignment*>::const_iterator m1 = map1.begin();
       m1 != map1.end(); ++m1) {
    keys1.push_back(m1->first);
  }
  list<int> keys2;
  for (map<int, NodeAlignment*>::const_iterator m2 = map2.begin();
       m2 != map2.end(); ++m2) {
    keys2.push_back(m2->first);
  }
  set_intersection(keys1.begin(), keys1.end(), keys2.begin(), keys2.end(),
		   back_inserter(*shared_keys));
}

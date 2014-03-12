#include <vector>
#include <string>
#include <sstream>

#include "src/AlignmentFilters.h"
#include "src/common.h"

using namespace std;

string GetCigarString(vector<BamTools::CigarOp>& cigar_ops){
  stringstream cigar_string;
  for(vector<BamTools::CigarOp>::iterator cigar_iter = cigar_ops.begin(); cigar_iter != cigar_ops.end(); cigar_iter++)
    cigar_string << cigar_iter->Type << cigar_iter->Length;
  return cigar_string.str();
}

template<typename CigarIterator> int GetDistToIndel(CigarIterator iter, CigarIterator end){
  if (iter != end && iter->Type == 'H')
    iter++;
  if (iter != end && iter->Type == 'S')
    iter++;
  if (iter != end){
    char type = iter->Type;
    if (type == 'M')
      return iter->Length;
    else if (type == 'I' || type == 'D')
      return 0;
    else if (type == 'S' || type == 'H')
      return 0;
    else 
      PrintMessageDieOnError("Invalid CIGAR char " + type, ERROR);
  }
  return 0;
}

pair<int,int> GetEndDistToIndel(AlignedRead& aln){
  vector<BamTools::CigarOp>::iterator cigar_iter;
  vector<BamTools::CigarOp>::iterator cigar_end;
  vector<int> vals;
  int head_dist = GetDistToIndel(aln.cigar_ops.begin(), aln.cigar_ops.end());
  int tail_dist = GetDistToIndel(aln.cigar_ops.rbegin(), aln.cigar_ops.rend());
  return pair<int,int>(head_dist, tail_dist);
}



pair<int,int> GetNumEndMatches(AlignedRead* aln, const string& ref_seq, int ref_seq_start){
  if (aln->read_start < ref_seq_start)
    return pair<int,int>(-1,0);

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
	return pair<int,int>(0, -1);
      if (read_index + cigar_iter->Length > aln->nucleotides.size())
	break;
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
    else 
      PrintMessageDieOnError("Invalid CIGAR char "+cigar_iter->Type, ERROR);
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
      return pair<int,int>(0,-1);
    else
      PrintMessageDieOnError("Improperly formatted CIGAR string", ERROR);
  }
  
  // Ensure that CIGAR string corresponded to aligned bases
  if (read_index != aln->nucleotides.size()){
    if (ref_index >= ref_seq.size())
      return pair<int,int>(0,-1);
    else
      PrintMessageDieOnError("CIGAR string does not correspond to alignment bases", ERROR);
  }

  if (beginning)
    return pair<int,int>(match_run, match_run);
  else
    return pair<int,int>(head_match, match_run);
}

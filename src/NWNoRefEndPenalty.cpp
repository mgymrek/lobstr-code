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

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

#include "common.h"				    
#include "NWNoRefEndPenalty.h"

namespace NWNoRefEndPenalty {

  // Create alignment scoring matrix
  const float a           =  2.0; // Match
  const float b           = -2.0; // Mismatch
  const float GAPOPEN     =  6.0;
  const float GAPEXTEND   =  0.125;
  const float s[4][4]     = {{ a, b, b, b },
			     { b, a, b, b },
			     { b, b, a, b },
			     { b, b, b, a }};

  // Large value used to penalize impossible configurations
  const float LARGE = 1000000;

  int base_to_int(char c){
    c = toupper(c);
    switch(c){
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    case 'N':
      return 4;
    default:
      PrintMessageDieOnError("Invalid character in read " + c, ERROR);
    }
    return -1;
  }


  float bestIndex(float s1, float s2, float s3, int* ptr){
    if (s2 > s1){
      if (s2 > s3){
	*ptr = 1;
	return s2;
      }
      else {
	*ptr = 2;
	return s3;
      }
    }
    else {
      if (s3 > s1){
	*ptr = 2;
	return s3;
      }
      else {
	*ptr = 0;
	return s1;
      }
    }
  }

  /* 
     Determines the optimal alignment ending point. Scans each of the 
     three scoring matrices for the highest-scoring alignment. Only considers
     the last row of each matrix.
  */
  void findOptimalStop(int& L1, int& L2, 
		       std::vector<float>& M, std::vector<float>& Iref, std::vector<float>& Iread, 
		       float& best_val, int& best_col, int& best_type){
    best_val   = -LARGE;
    best_col   = -1;
    best_type  = -1;
    int column = 0;
    for (int i = L2*(L1+1); i < (L2+1)*(L1+1); i++){
      if (M[i] >= best_val){
	best_val  = M[i];
	best_col  = column;
	best_type = 0;
      }
    
      if (Iref[i] > best_val){
	best_val  = Iref[i];
	best_col  = column;
	best_type = 1;
      }
    
      if (Iread[i] > best_val){
	best_val  = Iread[i];
	best_col  = column;
	best_type = 2;
      }
      column++;
    }
  }

  void nw_helper(std::vector<float>& M, std::vector<float>& Iref, std::vector<float>& Iread, 
		 std::vector<int>& traceM, std::vector<int>& traceIref, std::vector<int>& traceIread,
		 const std::string& refseq, const std::string& readseq, 
		 std::string& refseq_al, std::string& readseq_al, 
		 float* score, std::vector<BamTools::CigarOp>& cigar_list){
    int L1 = refseq.length();
    int L2 = readseq.length();
    cigar_list.clear();

    // Various variables used in the matrix calculations
    int ref_base, read_base, oindex, nindex;
    float s1, s2, s3;
    int c;

    // Fill in the 3 matrices using dynamic programming
    for (int i = 1; i <= L2; i++){
      for (int j = 1; j <= L1; j++){
	nindex    = i*(L1+1)+j;
	ref_base  = base_to_int(refseq[j-1]);
	read_base = base_to_int(readseq[i-1]);

	// Update M matrix (examine (i-1, j-1))
	oindex          = (i-1)*(L1+1)+(j-1);
	s1              = M[oindex];
	s2              = Iref[oindex];
	s3              = Iread[oindex];
	M[nindex]       = bestIndex(s1, s2, s3, &c) + s[ref_base][read_base];
	traceM[nindex]  = c;

	// Update Iref matrix (examine (i,j-1))
	oindex             = i*(L1+1) + (j-1);
	s1                 = M[oindex]     - GAPOPEN;
	s2                 = Iref[oindex]  - GAPEXTEND;
	s3                 = Iread[oindex] - GAPOPEN;
	Iref[nindex]       = bestIndex(s1, s2, s3, &c);
	traceIref[nindex]  = c;

	// Update Iread matrix (examine (i-1,j))
	oindex              = (i-1)*(L1+1) + j;
	s1                  = M[oindex]     - GAPOPEN;
	s2                  = Iref[oindex]  - GAPOPEN;
	s3                  = Iread[oindex] - GAPEXTEND;
	Iread[nindex]       = bestIndex(s1, s2, s3, &c);
	traceIread[nindex]  = c;
      }
    }
  
    //Find the best ending point for the alignment
    float best_val;
    int best_col, best_type;
    findOptimalStop(L1, L2, M, Iref, Iread, best_val, best_col, best_type);
  
    // Store the optimal alignment score
    *score = best_val;
  
    std::stringstream refseq_ss, readseq_ss, cigar_ss;
  
    // Handle trailing gaps
    for(int i = L1; i > best_col; i--){
      refseq_ss  << refseq.at(i-1);
      readseq_ss << "-";
    }

    // Traceback the optimal alignment
    int best_row = L2;
    std::string raw_cigar;
    int index;
    while (best_row > 0){
      index = best_row*(L1+1) + best_col;
      if (best_type == 0){
	// M
	refseq_ss  << refseq.at(best_col-1);
	readseq_ss << readseq.at(best_row-1);
	cigar_ss   << "M";
	best_type   = traceM[index];
	best_row--;
	best_col--;
      } 
      else if (best_type == 1){
	//Iref
	refseq_ss  << refseq.at(best_col-1);
	readseq_ss << "-";
	cigar_ss   << "D";
	best_type   = traceIref[index];
	best_col--;
      } 
      else if (best_type == 2){
	// Iread
	refseq_ss  << "-";
	readseq_ss << readseq.at(best_row-1);
	cigar_ss   << "I";
	best_type   = traceIread[index];
	best_row--;
      } 
      else
	PrintMessageDieOnError("Invalid matrix type in Needleman-Wunsch alignment", ERROR);
    }

    // Handle leading gaps
    for (int i = best_col; i > 0; i--){
      refseq_ss  << refseq.at(i-1);
      readseq_ss << "-";
    }
  
    // Order alignment front to back
    refseq_al  = refseq_ss.str();
    readseq_al = readseq_ss.str();
    raw_cigar  = cigar_ss.str();
    reverse(refseq_al.begin(),  refseq_al.end());
    reverse(readseq_al.begin(), readseq_al.end());
    reverse(raw_cigar.begin(),  raw_cigar.end());

    // Simplify cigar string
    char cigar_char = raw_cigar[0];
    int  num        = 1;
    char new_cigar_char;
    for(unsigned int i = 1; i < raw_cigar.length(); i++){
      new_cigar_char = raw_cigar[i];
      if (new_cigar_char != cigar_char){
	cigar_list.push_back(BamTools::CigarOp::CigarOp(cigar_char, num));
	num = 1;
	cigar_char = new_cigar_char;
      }
      else
	num += 1;
    }
    cigar_list.push_back(BamTools::CigarOp::CigarOp(cigar_char, num));
    if (cigar_list.back().Type == 'I')
      cigar_list.back().Type = 'S';
  }

  /* Trace back alignment using the string alingnments and the respective score matrices. */
  void traceAlignment(int L1, int L2, std::string ref_al, std::string read_al, 
		      std::vector<float>& M, std::vector<float>& Iref, std::vector<float>& Iread){
    std::cout << "Beginning trace..." 
	      << ref_al  << std::endl
	      << read_al << std::endl;
    int i = 0;
    int j = 0;
    for(unsigned int index = 0; index < ref_al.length(); index++){
      float score;
      if (ref_al[index] == '-'){
	i++;
	score = Iread[(L1+1)*i + j];
      }
      else if (read_al[index] == '-'){
	j++;
	score = Iref[(L1+1)*i + j];
      }
      else {
	i++;
	j++;
	score = M[(L1+1)*i + j];
      }
      std:: cout << index << " " << i << " " << j << " " << score << std::endl;
    }
    std::cout << "End of trace" << std::endl;
  }



  void Align(const std::string& ref_seq, const std::string& read_seq,
	     std::string& ref_seq_al, std::string& read_seq_al,
	     float* score, std::vector<BamTools::CigarOp>& cigar_list) {
    int L1       = ref_seq.length();
    int L2       = read_seq.length();
    int mat_size = (L1+1)*(L2+1);
  
    // Scoring matrices
    std::vector<float> M(mat_size);      // Ref and read bases aligned 
    std::vector<float> Iref(mat_size);   // Ref  base aligned with gap
    std::vector<float> Iread(mat_size);  // Read base aligned with gap

    // Traceback matrices
    std::vector<int> traceM(mat_size);
    std::vector<int> traceIref(mat_size);
    std::vector<int> traceIread(mat_size);

    M[0]     = 0.0;
    Iref[0]  = -LARGE; // Impossible
    Iread[0] = -LARGE; // Impossible

    // Fill in row (0,n)
    for(int i = 1; i < L1+1; i++){
      // No penalty for leading affine gap in reference seq 
      Iref[i]       = 0.0;  
      traceIref[i]  = 1;

      // Impossible
      Iread[i]      = -LARGE;
      traceIread[i] = -1;

      // Impossible
      M[i]      = -LARGE;
      traceM[i] = -1;
    }

    // Fill in column (n, 0)
    for(int i = 1; i < L2+1; i++){
      int index = i*(L1+1);

      // Penalty for leading affine gap in read sequence
      Iread[index]      = -GAPOPEN-(i-1)*GAPEXTEND;
      traceIread[index] = 2;

      // Impossible
      Iref[index]       = -LARGE; 
      traceIref[index]  = -1;

      // Impossible
      M[index]      = -LARGE;
      traceM[index] = -1;
    }
    
    // Determine optimal alignment
    nw_helper(M, Iref, Iread, traceM, traceIref, traceIread,
	      ref_seq, read_seq, ref_seq_al, read_seq_al, score, cigar_list);
  }
}

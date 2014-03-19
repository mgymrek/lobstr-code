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

#include <iostream>
#include <sstream>

#include <stdlib.h>

#include "src/common.h"
#include "src/ZAlgorithm.h"


namespace ZAlgorithm{
  void suffix_helper(const std::string& s1, const std::string& s2, int s2_left, int s2_right,
		     std::vector<int>& s1_matches, std::vector<int>& num_matches){
    num_matches  = std::vector<int>(s2_right - s2_left + 1, -1);
    int leftmost = s2_right+1, right_index = s2_right+1;
    for (int i = s2_right; i >= s2_left; i--){
      if (i <= leftmost){
	int index_a = s1.size()-1, index_b = i;
	while (index_a >= 0 && index_b >= 0 && s1[index_a] == s2[index_b]){
	  index_a--;
	  index_b--;
	}
	num_matches[i-s2_left] = i - index_b;
	if (index_b < i){
	  right_index = i;
	  leftmost    = index_b + 1;
	}
      }
      else {
	int twin     = i - right_index + s1.size()-1;
	int new_left = i - s1_matches[twin] + 1;
	if (new_left > leftmost)
	  num_matches[i-s2_left] = s1_matches[twin];
	else if (new_left < leftmost)
	  num_matches[i-s2_left] = i-leftmost+1;
	else {
	  int index_a = s1.size()-2-i+leftmost, index_b = leftmost-1;
	  while (index_a >= 0 && index_b >= 0 && s1[index_a] == s2[index_b]){
	    index_a--;
	    index_b--;
	  }
	  num_matches[i-s2_left] = i-index_b;
	  right_index            = i;
	  leftmost               = index_b + 1;
	}
      }
    }
  }

  void prefix_helper(const std::string& s1, const std::string& s2, int s2_left, int s2_right,
		     std::vector<int>& s1_matches, std::vector<int>& num_matches, int offset){
    num_matches = std::vector<int>(s2_right-s2_left+1+offset, -1);
    int rightmost = -1, left_index = -1;
    for (int i = s2_left; i <= s2_right; i++){
      if (i >= rightmost){
	int index_a = 0, index_b = i;
	while (index_a < static_cast<int>(s1.size()) && index_b < static_cast<int>(s2.size()) && s1[index_a] == s2[index_b]){
	  index_a++;
	  index_b++;
	}
	num_matches[i-s2_left+offset] = index_b - i;
	if (index_b > i){
	  left_index = i;
	  rightmost  = index_b - 1;
	}
      }
      else {
	int twin      = i - left_index;
	int new_right = i + s1_matches[twin] - 1;
	if (new_right < rightmost)
	  num_matches[i-s2_left+offset] = s1_matches[twin];
	else if (new_right > rightmost)
	  num_matches[i-s2_left+offset] = rightmost-i+1;
	else {
	  int index_a = rightmost+1-i, index_b = rightmost+1;
	  while (index_a < static_cast<int>(s1.size()) && index_b < static_cast<int>(s2.size()) && s1[index_a] == s2[index_b]){
	    index_a++;
	    index_b++;
	  }
	  num_matches[i-s2_left+offset] = index_b - i;
	  left_index     = i;
	  rightmost      = index_b - 1;
	}
      }
    }
  }


  void GetPrefixMatchCounts(const std::string& s1, const std::string& s2, std::vector<int>& num_matches) {
    std::vector<int> s1_matches;
    prefix_helper(s1, s1, 1, s1.size()-1, s1_matches, s1_matches,  1);
    prefix_helper(s1, s2, 0, s2.size()-1, s1_matches, num_matches, 0);
  }  

  void GetSuffixMatchCounts(const std::string& s1, const std::string& s2, std::vector<int>& num_matches) {
    std::vector<int> s1_matches;
    suffix_helper(s1, s1, 0, s1.size()-2, s1_matches, s1_matches);
    suffix_helper(s1, s2, 0, s2.size()-1, s1_matches, num_matches);
  }  

  void GetPrefixMatchCounts(const std::string& s1, const std::string& s2, int s2_start, int s2_stop, std::vector<int>& num_matches) {
    if (s2_start < 0 or s2_stop >= static_cast<int>(s2.size()))
	PrintMessageDieOnError("Invalid string indices provided to GetPrefixMatchCounts", ERROR);
    std::vector<int> s1_matches;
    prefix_helper(s1, s1, 1, s1.size()-1, s1_matches, s1_matches,  1);
    prefix_helper(s1, s2, s2_start, s2_stop, s1_matches, num_matches, 0);
  }  

  void GetSuffixMatchCounts(const std::string& s1, const std::string& s2, int s2_start, int s2_stop, std::vector<int>& num_matches) {
    if (s2_start < 0 or s2_stop >= static_cast<int>(s2.size()))
	PrintMessageDieOnError("Invalid string indices provided to GetSuffixMatchCounts", ERROR);
    std::vector<int> s1_matches;
    suffix_helper(s1, s1, 0, s1.size()-2, s1_matches, s1_matches);
    suffix_helper(s1, s2, s2_start, s2_stop, s1_matches, num_matches);
  }
}


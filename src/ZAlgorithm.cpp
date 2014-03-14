#include <iostream>
#include <sstream>

#include <stdlib.h>

#include "ZAlgorithm.h"


namespace ZAlgorithm{
  void suffix_helper(int start, string& s1, string& s2, 
		     vector<int>& s1_matches, vector<int>& num_matches){
    num_matches  = vector<int>(s2.size(), 0);
    int leftmost = s2.size(), right_index = s2.size();
    for (int i = start; i >= 0; i--){
      if (i <= leftmost){
	int index_a = s1.size()-1, index_b = i;
	while (index_a >= 0 && index_b >= 0 && s1[index_a] == s2[index_b]){
	  index_a--;
	  index_b--;
	}
	num_matches[i] = i - index_b;
	if (index_b < i){
	  right_index = i;
	  leftmost    = index_b + 1;
	}
      }
      else {
	int twin     = i - right_index + s1.size()-1;
	int new_left = i - s1_matches[twin] + 1;
	if (new_left > leftmost)
	  num_matches[i] = s1_matches[twin];
	else if (new_left < leftmost)
	  num_matches[i] = i-leftmost+1;
	else {
	  int index_a = s1.size()-2-i+leftmost, index_b = leftmost-1;
	  while (index_a >= 0 && index_b >= 0 && s1[index_a] == s2[index_b]){
	    index_a--;
	    index_b--;
	  }
	  num_matches[i] = i-index_b;
	  right_index    = i;
	  leftmost       = index_b + 1;
	}
      }
    }
  }


  void prefix_helper(unsigned int start, string& s1, string& s2, 
		     vector<int>& s1_matches, vector<int>& num_matches){
    num_matches = vector<int>(s2.size(), 0);
    int rightmost = 0, left_index = 0;
    for (int i = start; i < s2.size(); i++){
      if (i >= rightmost){
	int index_a = 0, index_b = i;
	while (index_a < s1.size() && index_b < s1.size() && s1[index_a] == s2[index_b]){
	  index_a++;
	  index_b++;
	}
	num_matches[i] = index_b - i;
	if (index_b > i){
	  left_index = i;
	  rightmost  = index_b - 1;
	}
      }
      else {
	int twin      = i - left_index;
	int new_right = i + s1_matches[twin] - 1;
	if (new_right < rightmost)
	  num_matches[i] = s1_matches[twin];
	else if (new_right > rightmost)
	  num_matches[i] = rightmost-i+1;
	else {
	  int index_a = rightmost+1-i, index_b = rightmost+1;
	  while (index_a < s1.size() && index_b < s2.size() && s1[index_a] == s2[index_b]){
	    index_a++;
	    index_b++;
	  }
	  num_matches[i] = index_b - i;
	  left_index     = i;
	  rightmost      = index_b - 1;
	}
      }
    }
  }

  void GetPrefixMatchCounts(string& s1, string& s2, vector<int>& num_matches) {
    vector<int> s1_matches;
    prefix_helper(1, s1, s1, s1_matches, s1_matches);
    prefix_helper(0, s1, s2, s1_matches, num_matches);
  }  

  void GetSuffixMatchCounts(string& s1, string& s2, vector<int>& num_matches) {
    vector<int> s1_matches;
    suffix_helper(s1.size()-2, s1, s1, s1_matches, s1_matches);
    suffix_helper(s2.size()-1, s1, s2, s1_matches, num_matches);
  }  
}


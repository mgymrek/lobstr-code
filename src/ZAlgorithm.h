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

#include <string>
#include <vector>


namespace ZAlgorithm{
  /* 
   * For each position in s2, calculates the length of the matching prefix of s1 and s2[i...]
   * and stores it in num_matches[i]. The provided vector is cleared and sized appropriately.
   * Runs in O(length_of_s1 + length_of_s2)
   */
  void GetPrefixMatchCounts(const std::string& s1, const std::string& s2, std::vector<int>& num_matches);
  
  /* 
   * For each position in s2, calculates the length of the matching suffix of s1 and s2[0...i]
   * and stores it in num_matches[i]. The provided vector is cleared and sized appropriately.
   * Runs in O(length_of_s1 + length_of_s2)
   */
  void GetSuffixMatchCounts(const std::string& s1, const std::string& s2, std::vector<int>& num_matches);

  /* 
   * For each position i in s2 in the range [s2_start, s2_stop], calculates the length of 
   * the matching prefix of s1 and s2[i...] and stores it in num_matches[i-s2_start]. 
   * The provided vector is cleared and sized appropriately.
   * Runs in O(length_of_s1 + size_of_s2_range)
   */
  void GetPrefixMatchCounts(const std::string& s1, const std::string& s2, int s2_start, int s2_stop, std::vector<int>& num_matches);
  
  /* 
   * For each position i in s2 in the range [s2_start, s2_stop], calculates the length of 
   * the matching suffix of s1 and s2[0...i] and stores it in num_matches[i-s2_start]. 
   * The provided vector is cleared and sized appropriately.
   * Runs in O(length_of_s1 + size_of_s2_range)
   */
  void GetSuffixMatchCounts(const std::string& s1, const std::string& s2, int s2_start, int s2_stop, std::vector<int>& num_matches);
}

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

using namespace std;


namespace ZAlgorithm{
  void GetPrefixMatchCounts(const string& s1, const string& s2, vector<int>& num_matches);
  void GetSuffixMatchCounts(const string& s1, const string& s2, vector<int>& num_matches);
}

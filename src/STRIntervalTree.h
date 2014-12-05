/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>

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

#ifndef SRC_STR_INTERVAL_TREE_H_
#define SRC_STR_INTERVAL_TREE_H_

#include <stdint.h>

#include <string>
#include <vector>

#include "src/IntervalTreeCore.h"
#include "src/ReferenceSTR.h"

struct Record {
  uint32_t start;
  uint32_t end;
};

using namespace std;

/*
Class to find intervals spanned by a given region
 */
class STRIntervalTree {
 public:
  STRIntervalTree();
  /* Load interval database */
  bool LoadIntervals(const std::vector<ReferenceSTR>& intervals);

  /* Get intervals spanned by a region */
  bool GetSpannedIntervals(const int& start, const int& end, std::vector<ReferenceSTR>* spanned_intervals);
  ~STRIntervalTree();

 private:
  IntervalTree<ReferenceSTR> itree;
};

#endif  // SRC_STR_INTERVAL_TREE_H_

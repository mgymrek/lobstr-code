/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>
Revisions     2014 Thomas Willems <twillems@mit.edu> 

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

// Much of this was shamelessly copied from BedTools
// using src/utils/BinTree/

#include <string.h>

#include "src/common.h"
#include "src/STRIntervalTree.h"

using namespace std;

STRIntervalTree::STRIntervalTree() {
}

bool STRIntervalTree::LoadIntervals(const vector<ReferenceSTR>& intervals) {
  vector<Interval<ReferenceSTR> > tree_intervals;
  for (vector<ReferenceSTR>::const_iterator it = intervals.begin();
       it != intervals.end(); it++) {
    tree_intervals.push_back(Interval<ReferenceSTR>(it->start, it->stop, *it));
  }
  itree = IntervalTree<ReferenceSTR>(tree_intervals);
  return true;
}

bool STRIntervalTree::GetSpannedIntervals(const int& start, const int& end, vector<ReferenceSTR>* spanned_intervals) {
  spanned_intervals->clear();
  vector<Interval<ReferenceSTR> > results;
  itree.findContained(start, end, results);
  for (vector<Interval<ReferenceSTR> >::iterator it = results.begin();
       it != results.end(); it++) {
    spanned_intervals->push_back(it->value);
  }
  return true;
}

bool STRIntervalTree::GetContainingRegions(const int& start, const int& end,
                                           vector<ReferenceSTR>* containing_intervals) {
  containing_intervals->clear();
  vector<Interval<ReferenceSTR> > results;
  itree.findOverlapping(start, end, results);
  for (vector<Interval<ReferenceSTR> >::iterator it=results.begin();
       it != results.end(); it++) {
    containing_intervals->push_back(it->value);
  }
  return true;
}

STRIntervalTree::~STRIntervalTree() {}

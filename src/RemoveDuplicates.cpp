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

#include <map>

#include "src/RemoveDuplicates.h"

using namespace std;

namespace RemoveDuplicates {

  void RemovePCRDuplicates(std::list<AlignedRead>* aligned_reads) {
    // map of <start pos, length> -> aligned reads list
    map<pair<int, size_t>, list<AlignedRead> > pcr_duplicates;
    // Group into duplicates
    for (list<AlignedRead>::const_iterator it = aligned_reads->begin();
	 it != aligned_reads->end(); it++) {
      pair<int, size_t> key(it->read_start, it->nucleotides.length());
      if (pcr_duplicates.find(key) != pcr_duplicates.end()) {
	pcr_duplicates.at(key).push_back(*it);
      } else {
	list<AlignedRead> pcr_dup_reads;
	pcr_dup_reads.push_back(*it);
	pcr_duplicates.insert(pair< pair<int, size_t>, list<AlignedRead> >
			      (key, pcr_dup_reads));
      }
    }
    // Choose one rep from each group
    list<AlignedRead> reads_after_rmdup;
    for (map<pair<int, size_t>, list<AlignedRead> >::const_iterator
	   it = pcr_duplicates.begin();
	 it != pcr_duplicates.end(); it++) {
      AlignedRead rep_read;
      GetRepRead(it->second, &rep_read);
      reads_after_rmdup.push_back(rep_read);
    }
    *aligned_reads = reads_after_rmdup;
  }

  void GetRepRead(const list<AlignedRead>& aligned_reads, AlignedRead* rep_alignment) {
    float highest_qual_score = 0;
    *rep_alignment = aligned_reads.front();
    for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
	 it != aligned_reads.end(); it++) {
      float qual = GetScore(it->qualities);
      if (qual > highest_qual_score) {
	*rep_alignment = *it;
      }
    }
  }

  float GetScore(const string& quality_string) {
    if (quality_string.empty()) return 0.0;
    float total_quality = 0;
    for (size_t i = 0; i < quality_string.length(); i++) {
      total_quality +=  quality_string.at(i) - 33;
    }
    return total_quality/static_cast<float>(quality_string.length());
  }
}

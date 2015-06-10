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

#include <err.h>

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "src/ReferenceSTR.h"
#include "src/runtime_parameters.h"

using namespace std;

int MAXCHUNKSPAN = 5000000;

ReferenceSTRContainer::ReferenceSTRContainer(const vector<ReferenceSTR>& _ref_strs) {
  ref_strs = _ref_strs;
  counter = 0;
  num_refs = ref_strs.size();
}

bool ReferenceSTRContainer::GetNextChunk(vector<ReferenceSTR>* ref_str_chunk,
                                         string* chrom, int* begin, int* end) {
  ref_str_chunk->clear();
  int current_chunk_size = 0;
  if (counter == num_refs) {
    return false;
  }
  *chrom = ref_strs.at(counter).chrom;
  *begin = ref_strs.at(counter).start;
  *end = ref_strs.at(counter).stop;
  while (current_chunk_size < CHUNKSIZE) {
    // Break if we're at the end of the list
    if (counter == num_refs) break;
    // don't span more than one chrom for a chunk
    if (ref_strs.at(counter).chrom != *chrom) break;
    // make sure this isn't going to make the chunk too big
    int mincoord = (ref_strs.at(counter).start < *begin)?ref_strs.at(counter).start:*begin;
    int maxcoord = (ref_strs.at(counter).stop > *end)?ref_strs.at(counter).stop:*end;
    if (maxcoord - mincoord > MAXCHUNKSPAN) break;
    // Else add to the list and update coords
    ref_str_chunk->push_back(ref_strs.at(counter));
    *begin = mincoord;
    *end = maxcoord;
    // Update counters
    current_chunk_size++;
    counter++;
  }
  return true;
}

bool ReferenceSTRContainer::GetChromChunk(vector<ReferenceSTR>* chrom_chunk,
					  const std::string& chrom) {
  size_t orig_size = chrom_chunk->size();
  for (vector<ReferenceSTR>::const_iterator it = ref_strs.begin();
       it != ref_strs.end(); it++) {
    if (it->chrom == chrom) {
      chrom_chunk->push_back(*it);
    }
  }
  size_t new_size = chrom_chunk->size();
  return (new_size > orig_size);
}


ReferenceSTRContainer::~ReferenceSTRContainer() {}

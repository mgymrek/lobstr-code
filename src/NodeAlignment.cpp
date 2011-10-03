/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <iostream>
#include <stdio.h>

#include "NodeAlignment.h"
#include "runtime_parameters.h"

using namespace std;

NodeAlignment::NodeAlignment() {
  // negative values mean it hasn't been initialized yet
  mismatches = -1;
}

NodeAlignment::~NodeAlignment() {}

void NodeAlignment::GetAlignmentWithID(int ident, AlignmentLocation *alnLoc){
  *alnLoc = string_id_to_alignment.at(ident);
}

void NodeAlignment::AddAlignmentLocation(const AlignmentLocation& alignment_location) {
  if (string_id_to_alignment.find(alignment_location.ident) !=
      string_id_to_alignment.end()) {
    // we already have an alignment to this string.
    // throw it away since it is a nonunique match
    string_id_to_alignment.erase(alignment_location.ident);
    string_identifiers_list.remove(alignment_location.ident);
    return;
  }
  string_identifiers_list.push_back(alignment_location.ident);
  string_id_to_alignment.insert(pair<int, AlignmentLocation>
				(alignment_location.ident,
				 alignment_location));
}

int NodeAlignment::GetStartLocation(const int& alignment_location_id) {
  return string_id_to_alignment.at(alignment_location_id).index;
}

int NodeAlignment::GetEndLocation(const int& alignment_location_id) {
  return string_id_to_alignment.at(alignment_location_id).index +
    len_aligned_string + 1;
}

void NodeAlignment::SetAlignedStringLength(int _string_len) {
  len_aligned_string = _string_len;
}

void NodeAlignment::GetStringIdentifiers(list<int>* identifiers) {
  *identifiers = string_identifiers_list;
}

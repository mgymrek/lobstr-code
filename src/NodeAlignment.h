/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <list>
#include <map>
#include <set>
using namespace std;

#ifndef NODEALIGNMENT_H_
#define NODEALIGNMENT_H_

// information about an alignment at a GST node
struct AlignmentLocation{
  // identity of string aligned to
  int ident;
  // index into string at this node
  int index;
};

class NodeAlignment {
 public:
  NodeAlignment();
  virtual ~NodeAlignment();

  // return the alignment location with for a given identifier
  void GetAlignmentWithID(int ident, AlignmentLocation *alnLoc);

  // add an alignment to a NodeAlignment object
  void AddAlignmentLocation(const AlignmentLocation& alignment_location);

  // get the start location of an alignment
  int GetStartLocation(const int& alignment_location_id);

  // get the end location of an alignment
  int GetEndLocation(const int& alignment_location_id);

  // set the length of the aligned string
  void SetAlignedStringLength(int _string_lewn);

  // map of sequence id to alignment locations
  // note discard if multiple alignments in the same string,
  // means we don't know where it came from
  map<int, AlignmentLocation> string_id_to_alignment;
  
  // get string identifiers
  void GetStringIdentifiers(list<int>* identifiers);

 private:
  // num mismatches in the alignment
  int mismatches;

  // length of the string at this node
  int len_aligned_string;
  
  // ids of the strings aligned to at this node
  list<int> string_identifiers_list;

};

#endif /* NODEALIGNMENT_H_ */

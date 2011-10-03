/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef READALIGNER_H_
#define READALIGNER_H_

#include <map>
#include <set>
#include <vector>

#include "GST.h"
#include "GSTContainer.h"
#include "MSReadRecord.h"

class ReadAligner {
 public:
  ReadAligner(map<string,GST*> *gstsL, map<string,GST*> *gstsR,
	      map<int,MSRecord>* msDict);
  virtual ~ReadAligner();
  // Main function: try to align the read to the suffix trees.
  // If we can align it, return true.
  bool ProcessRead(MSReadRecord* read);
 protected:
  // determine the STR repeat sequence from a string of nucleotides
  // get the canonicalized form of the repeat
  static bool getMSSeq(const string& nucs, int k, string* repeat);

  // Try aligning a read. Get  the L/R node alignments and the 
  // string identifier of an alignment.
  int AlignRead(const MSReadRecord& read, bool reverse,
		NodeAlignment* left_node_alignment,
		NodeAlignment* right_node_alignment,
		int* alignment_id, int* left_mismatch,
		int* right_mismatch);

  // From a list of all node alignments, get all identifiers that
  // occur only once. Keep track of which NodeAlignment is responsible
  // for which identifier
  static void GetUniqueIdentifiers(list<NodeAlignment>* node_alignments,
				   map<int, NodeAlignment*>* unique_identifiers);

  // From two maps, get keys that are shared between the two
  static void GetSharedKeys(const map<int, NodeAlignment*>& map1,
			    const map<int, NodeAlignment*>& map2,
			    list<int>* shared_identifiers);

  friend class ReadAlignerTest;

  map<string,GST*>* gstsL;
  map<string,GST*>* gstsR;
  map<int,MSRecord>* msDict;
};

#endif /* READALIGNER_H_ */

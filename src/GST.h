/*
 * GST.h
 *
 *  Created on: Apr 19, 2011
 *      Author: mgymrek
 *
 *      This class creates a generalized suffix tree based on a previous python implementation of suffix trees
 */

#ifndef GST_H_
#define GST_H_

#include <map>
#include <list>
#include <iostream>
#include "SuffixTree.h"
#include "MSRecord.h"
#include "GSTNode.h"
#include "NodeAlignment.h"
#include <string>


using namespace std;



class GST {
  public:
  GST(list<MSRecord> _msrecs, const int left);
  GST();
  virtual ~GST();

  // main function called externally that performs fuzzy matching
  list<NodeAlignment> fuzzyMatch(const string& str, int maxdist,
				 int maxalign, int numalign);  
  map<int,int> ids;

 private:
  // suffix tree holding the sequences
  SuffixTree* suffix_tree_;

  // internal function to call fuzzy match. fuzzy string matching
  // using DFS. match string starting from node,
  // maxdist away, return maximum maxalign alignments, numalign aligned
  // so far
  list<NodeAlignment> fuzzyMatch_int(const string& str,int maxdist,
				     SuffixNode* node, int maxalign,
				     int numalign);

  // check whether a node is a leaf
  bool isLeaf(SuffixNode* node);

  // get the hamming distance between a string and a char*
  static int hamming(const string& string1, const char* string2, int len2);

  // get the children nodes of a node
  list<SuffixNode*>  getChildren(SuffixNode* node);

  // walk the stree to get the node information of a single node
  // store in table so we don't have to calculate for a given node more than once
  NodeAlignment getNodeInfo(SuffixNode* node);

};

#endif /* GST_H_ */

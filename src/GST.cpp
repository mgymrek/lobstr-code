/*
 * GST.cpp
 *
 *  Created on: Apr 19, 2011
 *      Author: mgymrek
 */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "GST.h"
#include "runtime_parameters.h"

using namespace std;

const int kMaxReferenceSequenceLen = 1000;
const int kMinReferenceSequenceLen = 1;

GST::GST(list<MSRecord> _msrecs, const int left) {
  // if left == 1, then left flank, else right
  suffix_tree_ = new_SuffixTree(128); // TODO use DNA library?
  // count number of sequences added
  int counter = 0;
  // return code of adding string to suffix tree
  int addstr;
  for(list<MSRecord>::const_iterator it = _msrecs.begin(); it != _msrecs.end(); it++){
    if(it->leftFlank.size() > kMaxReferenceSequenceLen ||
       it->rightFlank.size() > kMaxReferenceSequenceLen ||
       it->leftFlank.size() < kMinReferenceSequenceLen ||
       it->rightFlank.size() < kMinReferenceSequenceLen){
      //      cerr << "Problem adding string to the tree... "
      //	   << " flanking regions too long or too short\n";
      //cerr << "left flank " << (*it).leftFlank.size() << " " <<
      //	"right flank " << (*it).rightFlank.size() << endl;
      continue;
    }
    // get sequence to add to suffix tree
    char* seq;
    if(left){seq = (char*)it->leftFlank.c_str();}
    else{seq = (char*)it->rightFlank.c_str();}
    addstr = SuffixTree_add(suffix_tree_,seq,counter);
    // check if adding was successful
    if (addstr == 0){
      cerr << "problem adding string to the tree... string " << seq << " was not added to the tree " << counter << endl;
    }else{
      ids.insert(pair<int,int>(counter,it->seqid));
    }
    counter++;
  }
}

GST::~GST() {
  ids.clear();
  delete_SuffixTree(suffix_tree_);
}

// external call to fuzzy matching
list<NodeAlignment> GST::fuzzyMatch(const string& str, int _maxdist,
				    int _maxalign, int _numalign){
  SuffixNode* root = SuffixTree_root(suffix_tree_);
  list<NodeAlignment> alignments;
  string str_to_align = str;
  alignments= fuzzyMatch_int(str, _maxdist, root, _maxalign, _numalign);
  free(root);
  return alignments;
}

// internal fuzzy matching
list<NodeAlignment> GST::fuzzyMatch_int(const string& str, int maxdist,
					SuffixNode* node, int maxalign,
					int numalign){
  // keep track of alignments
  list<NodeAlignment> alignments;

  // length of the string being aligned
  size_t string_to_align_len = str.size();

  // string leading to this node
  char* node_edgestr = SuffixNode_edgestr(node);

  // length of the edge string leading to this node
  int node_edgelen = SuffixNode_edgelen(node);

  // hamming distance between the string and the edgestring
  int node_hamming_dist = GST::hamming(str, node_edgestr, node_edgelen);

  // keep track of child alignments TODO move declaration later, right now for some
  // reason it explodes memory if I don't put it here
  list<SuffixNode*>children;
  list<NodeAlignment> newalns;
  int diff;
  // case 1: we're at the end of the string on this node
  if (string_to_align_len <= node_edgelen){
    if(node_hamming_dist <= maxdist){ // found alignment
      NodeAlignment node_alignment = getNodeInfo(node);
      node_alignment.SetAlignedStringLength(str.size());
      alignments.push_back(node_alignment);
    }
    free(node_edgestr);
    while(!children.empty()) free( children.front()),children.pop_front();
    /* if (align_debug) {
      cout << "case 1: end of string on this node" << endl;
      }*/
    return alignments;
  } 
  
  // case 2: we go past this node
  if(node_hamming_dist > maxdist){
    free(node_edgestr);
    while(!children.empty()) free(children.front()),children.pop_front();
    return alignments;
  } // didn't make it past this node

  // remaining part of string to be aligned in children nodes
  string newstr = str.substr(node_edgelen, string_to_align_len - node_edgelen);
  // length of new string to align
  int newstrlen = newstr.size();
  assert(newstrlen > 0);
  SuffixNode* first_child = SuffixNode_find_child(node, newstr[0]);
  int first_child_ident = -1;
  if (first_child != NULL) {
    first_child_ident = SuffixNode_ident(first_child);
    int num_remaining_mismatches = maxdist - node_hamming_dist;
    newalns = GST::fuzzyMatch_int(newstr, num_remaining_mismatches,
				  first_child, maxalign, numalign);
    if (newalns.size() > 0) {
      alignments.insert(alignments.end(),
			newalns.begin(),
			newalns.end());
      numalign = alignments.size();
      free(first_child);
      free(node_edgestr);
      return alignments;
    }
    free(first_child);
  }
  // we know there is at least one mismatch. If allowed mismatches is used up,
  // quit here.
  if (maxdist - node_hamming_dist < 0) {
    free(node_edgestr);
    /*	if (align_debug) {
	cout << "allowed mismatches used up" << endl;
	}*/
    return alignments;
  }
  
  // else it's on to the kids
  children = GST::getChildren(node);
  
  for(list<SuffixNode*>::iterator it = children.begin(); it != children.end(); it++){
    if (SuffixNode_ident(*it) == first_child_ident) continue;
    if (numalign <= maxalign){
      int num_remaining_mismatches = maxdist - node_hamming_dist;
      newalns = GST::fuzzyMatch_int(newstr, num_remaining_mismatches,
				    *it, maxalign, numalign);
      alignments.insert(alignments.end(),newalns.begin(),newalns.end());
      numalign += newalns.size();
    }
    newalns.clear(); // added
  }
  free(node_edgestr);
  while(!children.empty()) free (children.front()),children.pop_front();
  return alignments;
}

bool isLeaf(SuffixNode* node) {
  return (SuffixNode_num_leaves(node) == 1 && SuffixNode_num_children(node) == 0);
}

int GST::hamming(const string& string1, const char* string2, int len2){
  int hamming_distance = 0;
  int len_to_check = min((int)string1.size(),len2);
  for (int i = 0; i < len_to_check; ++i) {
    if ((string1[i] != string2[i])) {
      hamming_distance++;
    }
  }
  return hamming_distance;
}
list<SuffixNode*> GST::getChildren(SuffixNode* node){
  list<SuffixNode*> children;
  // if no children, you shall not pass
  int num_children = SuffixNode_num_children(node);
  if(num_children == 0){
    return children;
  }
  SuffixNode* current = SuffixNode_children(node);
  children.push_back(current);

  for(int i = 0; i < num_children - 1; i++){
    current = SuffixNode_next(current);
    children.push_back(current);
  }
  return children;
}

NodeAlignment GST::getNodeInfo(SuffixNode* node){
  NodeAlignment gst_node_alignment;
  // get any leaves at this node and update its table
  int num_leaves = SuffixNode_num_leaves(node);
  for(int i = 0; i < num_leaves; i++){
    // string ending at leaf
    char* tempstring;
    // index into string ending at leaf
    int index;
    // identity of leaf
    int ident;
    SuffixNode_leaf(node,i+1,&tempstring,&index,&ident);
    // set fields of alignment_location
    AlignmentLocation alignment_location;
    alignment_location.index = index;
    alignment_location.ident = ids.at(ident);
    gst_node_alignment.AddAlignmentLocation(alignment_location);
    free(tempstring);
  }

  // get children node information
  list<SuffixNode*> children = GST::getChildren(node);
  if (children.size() == 0) {return gst_node_alignment;}

  // get table for each of children and add it here
  int edgelen;
  for(list<SuffixNode*>::const_iterator it = children.begin();
      it != children.end(); it++){
    edgelen = SuffixNode_edgelen(*it);
    // get table for each of children and add it here
    NodeAlignment child_node_alignment;
      
    child_node_alignment = GST::getNodeInfo(*it);
    for(map<int, AlignmentLocation>::const_iterator it2 =
	  child_node_alignment.string_id_to_alignment.begin();
	it2 != child_node_alignment.string_id_to_alignment.end(); it2++){
      AlignmentLocation new_alignment_location;
      new_alignment_location.ident = (*it2).second.ident;
      new_alignment_location.index = (*it2).second.index;
      gst_node_alignment.AddAlignmentLocation(new_alignment_location);
    }
  }
  while(!children.empty()) free( children.front()),children.pop_front();
  return gst_node_alignment;
}



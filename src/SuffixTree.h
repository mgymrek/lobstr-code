#ifndef __SuffixTree_H__
#define __SuffixTree_H__

#include "stree.h"

/* A small utility class to make working with stree a little more
   pleasant. */

struct SuffixTree {
  STREE_STRUCT *tree;
};

struct SuffixNode {
  STREE_STRUCT *tree;
  SNODE_STRUCT *node;
};

typedef struct SuffixTree SuffixTree;
typedef struct SuffixNode SuffixNode;


/* SuffixTree methods */
SuffixTree* new_SuffixTree(int alphabet);
void delete_SuffixTree(SuffixTree* self);
int SuffixTree_add(SuffixTree *self, char* s, int id);
int SuffixTree_num_nodes(SuffixTree *self);
SuffixNode* SuffixTree_root(SuffixTree *self);
int SuffixTree_match(SuffixTree *self, char *s,
		     SuffixNode **node_out, int *endpos_out);


/* SuffixNode methods */
void delete_SuffixNode(SuffixNode *self);
int SuffixNode_num_children(SuffixNode *self);
SuffixNode* SuffixNode_find_child(SuffixNode *self, char ch);
SuffixNode* SuffixNode_children(SuffixNode *self);
SuffixNode* SuffixNode_next(SuffixNode *self);
SuffixNode* SuffixNode_parent(SuffixNode *self);
SuffixNode* SuffixNode_suffix_link(SuffixNode *self);
int SuffixNode_edgelen(SuffixNode *self);
char* SuffixNode_edgestr(SuffixNode *self);
char* SuffixNode_getch(SuffixNode *self);
int SuffixNode_labellen(SuffixNode *self);
char* SuffixNode_labelstr(SuffixNode *self);
int SuffixNode_ident(SuffixNode *self);
int SuffixNode_num_leaves(SuffixNode *self);
int SuffixNode_leaf(SuffixNode *self, int leafnum,
		    char **string_out, int *pos_out, int *id_out);



#endif

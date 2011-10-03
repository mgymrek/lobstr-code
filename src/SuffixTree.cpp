#include "SuffixTree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
/**********************************************************************/

SuffixTree* new_SuffixTree(int alphabet) {
  const int COPY_FLAG = 1;
  SuffixTree *self = (SuffixTree *)malloc(sizeof(SuffixTree));
  self->tree = stree_new_tree(alphabet, COPY_FLAG);
  return self;
}

void delete_SuffixTree(SuffixTree* self) {
  stree_delete_tree(self->tree);
  free(self);
}

int SuffixTree_add(SuffixTree *self, char* s, int id) {
  return stree_add_string(self->tree, s, strlen(s), id);
}


int SuffixTree_match(SuffixTree *self, char *s,
		     SuffixNode **node_out, int *endpos_out) {
  int matched_chars;
  STREE_NODE stree_node;
  matched_chars = stree_match(self->tree, s, strlen(s),
			      &stree_node, endpos_out);
  (*node_out) = (SuffixNode *)malloc(sizeof(SuffixNode *));
  (*node_out)->tree = self->tree;
  (*node_out)->node = stree_node;
  if ((*node_out)->node == NULL) {
    free(*node_out);
    (*node_out) = NULL;
  }
  return matched_chars;
}


int SuffixTree_num_nodes(SuffixTree *self) {
  return stree_get_num_nodes(self->tree);
}

SuffixNode* SuffixTree_root(SuffixTree *self) {
  SuffixNode *s = (SuffixNode *)malloc(sizeof(SuffixNode));
  s->tree = self->tree;
  s->node = stree_get_root(self->tree);
  if (s->node == NULL) {
    free(s);
    return NULL;
  }
  return s;
}



/**********************************************************************/

void delete_SuffixNode(SuffixNode *self) {
  free(self);
}


int SuffixNode_num_children(SuffixNode *self) {
  return stree_get_num_children(self->tree, self->node);
}


SuffixNode* SuffixNode_find_child(SuffixNode *self, char ch) {
  SuffixNode* child = (SuffixNode *)malloc(sizeof(SuffixNode));
  child->tree = self->tree;
  child->node = stree_find_child(self->tree, self->node, ch);
  if (child->node == NULL) {
    free(child);
    return NULL;
  }
  return child;
}


// We may want to change this so that it returns a list of children, so that
// we don't have to think about linked list stuff.
SuffixNode* SuffixNode_children(SuffixNode *self) {
  SuffixNode* child = (SuffixNode *)malloc(sizeof(SuffixNode));
  child->tree = self->tree;
  child->node = stree_get_children(self->tree, self->node);
  if (child->node == NULL) {
    free(child);
    return NULL;
  }
    return child;
}


SuffixNode* SuffixNode_next(SuffixNode *self) {
  SuffixNode* sibling = (SuffixNode *)malloc(sizeof(SuffixNode));
  sibling->tree = self->tree;
  sibling->node = stree_get_next(self->tree, self->node);
  if (sibling->node == NULL) {
    free(sibling);
    return NULL;
  }
  return sibling;
}


SuffixNode* SuffixNode_parent(SuffixNode *self) {
  SuffixNode* parent = (SuffixNode *)malloc(sizeof(SuffixNode));
  parent->tree = self->tree;
  parent->node = stree_get_parent(self->tree, self->node);
  if (parent->node == NULL) {
    free(parent);
    return NULL;
  }
  return parent;
}


SuffixNode* SuffixNode_suffix_link(SuffixNode *self) {
  SuffixNode* suffix = (SuffixNode *)malloc(sizeof(SuffixNode));
  suffix->tree = self->tree;
  suffix->node = stree_get_suffix_link(self->tree, self->node);
  if (suffix->node == NULL) {
    free(suffix);
    return NULL;
  }
  return suffix;
}


int SuffixNode_edgelen(SuffixNode *self) {
  return stree_get_edgelen(self->tree, self->node);
}


char* SuffixNode_edgestr(SuffixNode *self) {
  char *result_str = NULL;
  int len = stree_get_edgelen(self->tree, self->node);
  result_str = (char*)malloc(sizeof(char) * (len+1));
  strncpy(result_str, stree_get_edgestr(self->tree, self->node), len);
  result_str[len] = '\0';
  return result_str;
};


char* SuffixNode_getch(SuffixNode *self) {
  /** If the node is a root, it has no 'ch', so we must take direct
      action to avoid a segfault.
      SWIG.**/
  char *result;
  if (self->node == stree_get_root(self->tree)) {
    return NULL;
  }
  result = (char*)malloc(sizeof(char) + 1);
  result[0] = stree_getch(self->tree, self->node); result[1] = '\0';
  return result;
};


int SuffixNode_labellen(SuffixNode *self) {
  return stree_get_labellen(self->tree, self->node);
}


char* SuffixNode_labelstr(SuffixNode *self) {
  char *result_str = NULL;
  int len = stree_get_labellen(self->tree, self->node);
  result_str = (char*)malloc(sizeof(char) * (len+1));
  stree_get_label(self->tree, self->node, result_str,
		  len, 0);
  result_str[len] = '\0';   // Perhaps this isn't necessary...
  return result_str;
}


int SuffixNode_ident(SuffixNode *self) {
  return stree_get_ident(self->tree, self->node);
}


int SuffixNode_num_leaves(SuffixNode *self) {
  return stree_get_num_leaves(self->tree, self->node);
}


int SuffixNode_leaf(SuffixNode *self, int leafnum,
		    char **string_out, int *pos_out, int *id_out) {
  int result;
  char *temp_string;
  int len;
  result = stree_get_leaf(self->tree, self->node, leafnum,
			  &temp_string, pos_out, id_out);
  len = strlen(temp_string);
  (*string_out) = (char*)malloc(sizeof(char) * (len + 1));
  strncpy((*string_out), temp_string, len);
  (*string_out)[len] = '\0';
  return result;
}

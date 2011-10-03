/*
 * stree.c
 *
 * Implementation of the suffix tree data structure.
 *
 * NOTES:
 *    7/94  -  Initial implementation of Weiner's algorithm (Joyce Lau)
 *    9/94  -  Partial Implementation.  (James Knight)
 *    9/94  -  Completed the implementation.  (Sean Davis)
 *    9/94  -  Redid Sean's implementation.  (James Knight)
 *    9/95  -  Reimplemented suffix trees, optimized the data structure,
 *             created streeopt.[ch]   (James Knight)
 *    4/96  -  Modularized the code  (James Knight)
 *    2/00  -  Fixed a bug in function free_element() reported by Gregory R.
 *             Parker: increment page->count only once, not twice (Jens Stoye)
 *    2/01  -  Fixed a bug in function int_stree_set_idents() reported by
 *             Haidong Wang: clean the "dirty" flag (Jens Stoye)
 *    2/01  -  Fixed another bug in function free_element() reported by Haidong
 *             Wang: set pointer back->next correctly (Jens Stoye)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "stree.h"
using namespace std;
#include <iostream>

static int initmaps = 0;
static char dnamap[128], proteinmap[128], selfmap[128];

static void int_stree_init_alphamaps(void);

static int max_alloc = 0;
int cur_alloc = 0;

/*
 *
 * The Suffix Tree Interface Procedures
 *
 *
 */

/*
 * stree_new_tree
 *
 * Allocates a new suffix tree data structure.
 *
 * Parameters:  alphabet   -  the size of the alphabet
 *                               (can either be STREE_DNA, STREE_RNA,
 *                                STREE_PROTEIN or the size of the alphabet.
 *                                In the last case, all strings are assumed
 *                                to use "characters" in the range of
 *                                0..alphasize-1.)
 *              copyflag   -  make a copy of each sequence?
 *
 * Returns:  A SUFFIX_TREE structure
 */
SUFFIX_TREE stree_new_tree(int alphabet, int copyflag)
{
  SUFFIX_TREE tree;

  if (alphabet > MAXALPHA ||
      (alphabet <= 0 && alphabet != STREE_DNA && alphabet != STREE_PROTEIN))
    return NULL;

  if (!initmaps)
    int_stree_init_alphamaps();

  /*
   * Allocate the space.
   */
  if ((tree = (STREE_STRUCT*)malloc(sizeof(STREE_STRUCT))) == NULL)
    return NULL;
  memset(tree, 0, sizeof(STREE_STRUCT));

  if ((tree->strings = (char**)malloc(MAXNUMSTR * sizeof(char *))) == NULL ||
      (tree->lengths = (int*)malloc(MAXNUMSTR * sizeof(int))) == NULL ||
      (tree->ids = (int*)malloc(MAXNUMSTR * sizeof(int))) == NULL ||
      (tree->root = int_stree_new_node(tree, NULL, 0)) == NULL) {
    if (tree->strings != NULL)  free(tree->strings);
    if (tree->lengths != NULL)  free(tree->lengths);
    if (tree->ids != NULL)  free(tree->ids);
    if (tree->root != NULL) int_stree_free_node(tree, tree->root);
    free(tree);
    return NULL;
  }
  memset(tree->strings, 0, MAXNUMSTR * sizeof(char *));

  tree->num_nodes = 1;
  tree->copyflag = copyflag;
  if (alphabet == STREE_DNA) {
    tree->alpha_size = 4;
    tree->alpha_map = dnamap;
  }
  else if (alphabet == STREE_PROTEIN) {
    tree->alpha_size = 20;
    tree->alpha_map = proteinmap;
  }
  else {
    tree->alpha_size = alphabet;
    tree->alpha_map = selfmap;
  }

  tree->root->nextisparent = 1;

  return tree;
}


/*
 * stree_delete_tree
 *
 * Frees the SUFFIX_TREE data structure and all of its allocated space.
 *
 * Parameters:  tree  -  a suffix tree
 *
 * Returns:  nothing.
 */
void stree_delete_tree(SUFFIX_TREE tree)
{
  int i;

  int_stree_delete_subtree(tree, stree_get_root(tree));

  if (tree->strings != NULL) {
    if (tree->copyflag) {
      for (i=0; i < MAXNUMSTR; i++)
        if (tree->strings[i] != NULL)
          free(tree->strings[i]);
    }
    free(tree->strings);
  }
  if (tree->ids != NULL)
    free(tree->ids);
  if (tree->lengths != NULL)
    free(tree->lengths);

  free(tree);
}


/*
 * stree_add_string
 *
 * Implements Ukkonen's construction algorithm to add a string to the
 * suffix tree.
 *
 * This operation is an "atomic" operation.  In the far too likely case
 * that the program runs out of memory (or hits the maximum allocated
 * memory set by stree_set_max_alloc), this operation undoes any partial
 * changes it may have made to the algorithm, and it leaves the tree in
 * its original form.  Thus, you can just keep adding strings until the
 * function returns 0, and not have too worry about whether a call to the
 * function will trash the tree just because there's no memory left.
 *
 * NOTE:  The `id' value given must be unique to any of the strings
 *        added to the tree, and must be a small integer greater than
 *        0.
 *
 *        The best id's to use are to number the strings from 1 to K.
 *
 * Parameters:  tree  -  a suffix tree
 *              S     -  the sequence
 *              M     -  the sequence length
 *              id    -  the sequence identifier
 *              Sraw  -  the raw sequence (i.e. whose characters
 *                       are not translated to 0..alphasize-1)
 *
 * Returns:  non-zero on success, zero on error.
 */
int stree_add_string(SUFFIX_TREE tree, char *S, int M, int strid)
{

  int i, j, g, h, gprime, edgelen, id;
  char *edgestr;
  STREE_NODE node, lastnode, root, child, parent;
  STREE_LEAF leaf;
  id = int_stree_insert_string(tree, S, M, strid);
  if (id == -1)
    return 0;

  /*
   * Run Ukkonen's algorithm to add the string to the suffix tree.
   */
  root = stree_get_root(tree);
  node = lastnode = root;
  g = 0;
  edgelen = 0;
  edgestr = NULL;

  for (i=0,j=0; i <= M; i++)  {
    for ( ; j <= i && j < M; j++) {
      /*
       * Perform the extension from S[j..i-1] to S[j..i].  One of the
       * following two cases holds:
       *    a) g == 0, node == root and i == j.
       *         (meaning that in the previous outer loop,
       *          all of the extensions S[1..i-1], S[2..i-1], ...,
       *          S[i-1..i-1] were done.)
       *    b) g > 0, node != root and the string S[j..i-1]
       *       ends at the g'th character of node's edge.
       */
      if (g == 0 || g == edgelen) {
        if (i < M) {
          if ((child = stree_find_child(tree, node, S[i])) != NULL) {
            node = child;
            g = 1;
            edgestr = stree_get_edgestr(tree, node);
            edgelen = stree_get_edgelen(tree, node);
            break;
          }

          if ((leaf = int_stree_new_leaf(tree, id, i)) == NULL ||
              tree->num_nodes == MAXNUMNODES ||
              (node = int_stree_connect(tree, node,
                                        (STREE_NODE) leaf)) == NULL) {
            if (leaf != NULL)
              int_stree_free_leaf(tree, leaf);
            int_stree_remove_to_position(tree, id, j);
            int_stree_delete_string(tree, id);
            return 0;
          }

          tree->num_nodes++;
        }
        else {
          if ((int_stree_isaleaf(tree, node) &&
               (node = int_stree_convert_leafnode(tree, node)) == NULL)) {
            int_stree_remove_to_position(tree, id, j);
            int_stree_delete_string(tree, id);
            return 0;
          }

          if (!int_stree_add_intleaf(tree, node, id, j)) {
            int_stree_remove_to_position(tree, id, j);
            int_stree_delete_string(tree, id);
            return 0;
          }
        }

        if (lastnode != root && lastnode->suffix_link == NULL)
          lastnode->suffix_link = node;
        lastnode = node;
      }
      else {
        /*
         * g > 0 && g < edgelen, and so S[j..i-1] ends in the middle
         * of some edge.
         *
         * If the next character in the edge label matches the next
         * input character, keep moving down that edge.  Otherwise,
         * split the edge at that point and add a new leaf for the
         * suffix.
         */
        if (i < M &&
            stree_mapch(tree, S[i]) == stree_mapch(tree, edgestr[g])) {
          g++;
          break;
        }

        if (tree->num_nodes == MAXNUMNODES ||
            (node = int_stree_edge_split(tree, node, g)) == NULL) {
          int_stree_remove_to_position(tree, id, j);
          int_stree_delete_string(tree, id);
          return 0;
        }

        edgestr = stree_get_edgestr(tree, node);
        edgelen = stree_get_edgelen(tree, node);

        if (i < M) {
          if ((leaf = int_stree_new_leaf(tree, id, i)) == NULL ||
              tree->num_nodes == MAXNUMNODES ||
              (node = int_stree_connect(tree, node,
                                        (STREE_NODE) leaf)) == NULL) {
            if (leaf != NULL)
              int_stree_free_leaf(tree, leaf);
            int_stree_remove_to_position(tree, id, j);
            int_stree_delete_string(tree, id);
            return 0;
          }

          tree->num_nodes++;
        }
        else {
          if ((int_stree_isaleaf(tree, node) &&
               (node = int_stree_convert_leafnode(tree, node)) == NULL)) {
            int_stree_remove_to_position(tree, id, j);
            int_stree_delete_string(tree, id);
            return 0;
          }

          if (!int_stree_add_intleaf(tree, node, id, j)) {
            int_stree_remove_to_position(tree, id, j);
            int_stree_delete_string(tree, id);
            return 0;
          }
        }

        if (lastnode != root && lastnode->suffix_link == NULL)
          lastnode->suffix_link = node;
        lastnode = node;
      }

      /*
       * Now, having extended S[j..i-1] to S[j..i] by rule 2, find where
       * S[j+1..i-1] is.
       */
      if (node == root)
        ;
      else if (g == edgelen && node->suffix_link != NULL) {
        node = node->suffix_link;
        edgestr = stree_get_edgestr(tree, node);
        edgelen = stree_get_edgelen(tree, node);
        g = edgelen;
      }
      else {
        parent = stree_get_parent(tree, node);
        if (parent != root)
          node = parent->suffix_link;
        else {
          node = root;
          g--;
        }
        edgelen = stree_get_edgelen(tree, node);

        h = i - g;
        while (g > 0) {
          node = stree_find_child(tree, node, S[h]);
          gprime = stree_get_edgelen(tree, node);
          if (gprime > g)
            break;

          g -= gprime;
          h += gprime;
        }

        edgelen = stree_get_edgelen(tree, node);
        edgestr = stree_get_edgestr(tree, node);

        if (g == 0) {
          if (lastnode != root && !int_stree_isaleaf(tree, node) &&
              lastnode->suffix_link == NULL) {
            lastnode->suffix_link = node;
            lastnode = node;
          }

          if (node != root)
            g = edgelen;
        }
      }
    }
  }

  return 1;
}


/*
 * stree_remove_string
 *
 * Removes a string from the suffix tree, pruning branches and suffix
 * links from the tree and compacting the tree as necessary.
 *
 * Parameters:  tree   -  A suffix tree
 *              strid  -  The identifier of the sequence to be removed.
 *
 * Returns:  non-zero on success, zero on error.
 */
int stree_remove_string(SUFFIX_TREE tree, int strid)
{
  int i;

  for (i=0; i < MAXNUMSTR; i++) {
    if (tree->strings[i] != NULL && tree->ids[i] == strid)
      break;
  }

  if (i == MAXNUMSTR)
    return 0;

  int_stree_remove_to_position(tree, i, tree->lengths[i]);
  /* debug --- this is causing segfaults!  None of the strings are
     guaranteed to be null terminated, so maybe there's some problems
     here.

     I'm not convinced int_stree_remove_to_position is doing the
     right thing, as the code isn't used anywhere else except when
     stree_add_string fails... and that probably means that this code has
     not been tested!
  */
  int_stree_delete_string(tree, i);
  return 1;
}


/*
 * stree_traverse & stree_traverse_subtree
 *
 * Use a non-recursive traversal of the tree (or a subtree), calling the
 * two function parameters before and after recursing at each node, resp.
 * When memory is at a premium, this traversal may be useful.
 *
 * Note that either of the function parameters can be NULL, if you just
 * need to do pre-order or post-order processing.
 *
 * The traversal uses the `ch' field of the tree nodes to hold its
 * state information.  After the traversal is finished with a node, it
 * will restore that ch value.
 *
 * Parameters:  tree          -  a suffix tree
 *              node          -  root node of the traversal
 *                                 (stree_traverse_subtree only)
 *              preorder_fn   -  function to call before visiting the children
 *              postorder_fn  -  function to call after visiting all children
 *
 * Returns:  nothing.
 */

/*
 * stree_match & stree_walk
 *
 * Traverse the path down the tree whose path label matches T, and return
 * the number of characters of T matches, and the node and position along
 * the node's edge where the matching to T ends.
 *
 * Parameters:  tree      -  a suffix tree
 *              node      -  what node to start the walk down the tree
 *              pos       -  position along node's edge to start the walk
 *                              (`node' and `pos' are stree_walk only)
 *              T         -  the sequence to match
 *              N         -  the sequence length
 *              node_out  -  address of where to store the node where
 *                           the traversal ends
 *              pos_out   -  address of where to store the character position
 *                           along the ending node's edge of the endpoint of
 *                           the traversal
 *
 * Returns:  The number of characters of T matched.
 */
int stree_match(SUFFIX_TREE tree, char *T, int N,
                STREE_NODE *node_out, int *pos_out)
{
  return stree_walk(tree, stree_get_root(tree), 0, T, N, node_out, pos_out);
}

int stree_walk(SUFFIX_TREE tree, STREE_NODE node, int pos, char *T, int N,
               STREE_NODE *node_out, int *pos_out)
{
  int len, endpos, edgelen;
  char *edgestr;
  STREE_NODE endnode;

  len = int_stree_walk_to_leaf(tree, node, pos, T, N, &endnode, &endpos);

  if (!int_stree_isaleaf(tree, endnode) || len == N) {
    *node_out = endnode;
    *pos_out = endpos;
    return len;
  }

  edgestr = stree_get_edgestr(tree, endnode);
  edgelen = stree_get_edgelen(tree, endnode);

  while (len < N && endpos < edgelen &&
         stree_mapch(tree, T[len]) == stree_mapch(tree, edgestr[endpos])) {
    len++;
    endpos++;
  }

  *node_out = endnode;
  *pos_out = endpos;
  return len;
}


/*
 * stree_find_child
 *
 * Find the child of a node whose edge label begins with the character given
 * as a parameter.
 *
 * Parameters:  tree  -  a suffix tree
 *              node  -  a tree node
 *              ch    -  a character
 *
 * Returns:  a tree node or NULL.
 */
STREE_NODE stree_find_child(SUFFIX_TREE tree, STREE_NODE node, char ch)
{
  char mapch;
  STREE_NODE child;

  if (ch < 0 || ch >= tree->alpha_size)
    return NULL;

  mapch = stree_mapch(tree, ch);

  child = stree_get_children(tree, node);
  while (child != NULL && stree_get_mapch(tree, child) < mapch)
    child = stree_get_next(tree, child);

  if (child != NULL && mapch == stree_get_mapch(tree, child))
    return child;
  else
    return NULL;
}


/*
 * stree_get_num_children
 *
 * Return the number of children of a node.
 *
 * Parameters:  tree  -  a suffix tree
 *              node  -  a tree node
 *
 * Returns:  the number of children.
 */
int stree_get_num_children(SUFFIX_TREE tree, STREE_NODE node)
{
  int count;
  STREE_NODE child;

  count = 0;
  child = stree_get_children(tree, node);
  while (child != NULL) {
    count++;
    child = stree_get_next(tree, child);
  }

  return count;
}


/*
 * stree_get_labellen
 *
 * Get the length of the string labelling the path from the root to
 * a tree node.
 *
 * Parameters:  tree  -  a suffix tree
 *              node  -  a tree node
 *
 * Returns:  the length of the node's label.
 */
int stree_get_labellen(SUFFIX_TREE tree, STREE_NODE node)
{
  int len;

  len = 0;
  while (node != stree_get_root(tree)) {
    len += stree_get_edgelen(tree, node);
    node = stree_get_parent(tree, node);
  }
  return len;
}


/*
 * stree_get_label
 *
 * Get the string labelling the path from the root to a tree node and
 * store that string (or a part of the string) in the given buffer.
 *
 * If the node's label is longer than the buffer, then `buflen'
 * characters from either the beginning or end of the label (depending
 * on the value of `endflag') are copied into the buffer and the string
 * is NOT NULL-terminated.  Otherwise, the string will be NULL-terminated.
 *
 * Parameters:  tree     -  a suffix tree
 *              node     -  a tree node
 *              buffer   -  the character buffer
 *              buflen   -  the buffer length
 *              endflag  -  copy from the end of the label?
 *
 * Returns:  nothing.
 */
void stree_get_label(SUFFIX_TREE tree, STREE_NODE node, char *buffer,
                     int buflen, int endflag)
{
  int len, skip, edgelen;
  char *edgestr, *bufptr;

  len = stree_get_labellen(tree, node);
  skip = 0;

  if (buflen > len)
    buffer[len] = '\0';
  else {
    if (len > buflen && !endflag)
      skip = len - buflen;
    len = buflen;
  }

  /*
   * Fill in the buffer from the end to the beginning, as we move up
   * the tree.  If `endflag' is false and the buffer is smaller than
   * the label, then skip past the "last" `len - buflen' characters (i.e., the
   * last characters on the path to the node, but the first characters
   * that will be seen moving up to the root).
   */
  bufptr = buffer + len;
  while (len > 0 && node != stree_get_root(tree)) {
    edgelen = stree_get_edgelen(tree, node);

    if (skip >= edgelen)
      skip -= edgelen;
    else {
      if (skip > 0) {
        edgelen -= skip;
        skip = 0;
      }
      edgestr = stree_get_edgestr(tree, node) + edgelen;
      for ( ; len > 0 && edgelen > 0; edgelen--,len--)
        *--bufptr = *--edgestr;
    }

    node = stree_get_parent(tree, node);
  }
}


/*
 * stree_get_num_leaves
 *
 * Return the number of suffices that end at the tree node (these are the
 * "leaves" in the theoretical suffix tree).
 *
 * Parameters:  tree  -  a suffix tree
 *              node  -  a tree node
 *
 * Returns:  the number of "leaves" at that node.
 */
int stree_get_num_leaves(SUFFIX_TREE tree, STREE_NODE node)
{
  int count;
  STREE_INTLEAF intleaf;

  if (int_stree_isaleaf(tree, node))
    return 1;

  count = 0;
  intleaf = int_stree_get_intleaves(tree, node);
  while (intleaf != NULL) {
    count++;
    intleaf = intleaf->next;
  }

  return count;
}


/*
 * stree_get_leaf
 *
 * Get the sequence information about one of the suffices that end at
 * a tree node.  The `leafnum' parameter gives a number between 1 and
 * the number of "leaves" at the node, and information about that "leaf"
 * is returned.
 *
 * NOTE:  There is no ordering of the "leaves" at a node, either by
 *        string identifier number or by position.  You get whatever
 *        order the construction algorithm adds them to the node.
 *
 * Parameters:  tree        -  a suffix tree
 *              node        -  a tree node
 *              leafnum     -  which leaf to return
 *              string_out  -  address where to store the suffix pointer
 *                                (the value stored there points to the
 *                                 beginning of the sequence containing the
 *                                 suffix, not the beginning of the suffix)
 *              pos_out     -  address where to store the position of the
 *                             suffix in the sequence
 *                                (a value between 1 and the sequence length)
 *              id_out      -  address where to store the seq. identifier
 *
 * Returns:  non-zero if a leaf was returned (i.e., the `leafnum' value
 *           referred to a valid leaf), and zero otherwise.
 *           NOTE: If no leaf is returned, *string_out, *pos_out and *id_out
 *                 are left untouched.
 */
int stree_get_leaf(SUFFIX_TREE tree, STREE_NODE node, int leafnum,
                   char **string_out, int *pos_out, int *id_out)
{
  int i;
  STREE_INTLEAF intleaf;
  STREE_LEAF leaf;

  if (int_stree_isaleaf(tree, node)) {
    if (leafnum != 1)
      return 0;

    leaf = (STREE_LEAF) node;
    if (string_out) *string_out = tree->strings[leaf->strid];
    if (pos_out) *pos_out = int_stree_get_leafpos(tree, leaf);
    if (id_out) *id_out = tree->ids[leaf->strid];
    return 1;
  }
  else {
    intleaf = int_stree_get_intleaves(tree, node);
    for (i=1; intleaf != NULL && i < leafnum; i++,intleaf=intleaf->next) ;

    if (intleaf == NULL)
      return 0;
    else {
      if (string_out) *string_out = tree->strings[intleaf->strid];
      if (pos_out) *pos_out = intleaf->pos;
      if (id_out) *id_out = tree->ids[intleaf->strid];
      return 1;
    }
  }
}


/*
 * stree_set_max_alloc
 *
 * Sets the maximum number of bytes to use while constructing the
 * suffix trees.  This is a global count across all suffix trees being
 * constructed.
 *
 * Calling this function with a value of 0 removes any maximum bound
 * on the allocation, and the functions will allocate until they run
 * out of memory.  (This is the default.)
 *
 * NOTE:  This count does NOT include the space for the SUFFIX_TREE
 *        structures or the copied sequences.  It only counts the
 *        space used by the suffix trees themselves.
 *
 * Parameters:  size  -  the max allocation size or 0.
 *
 * Returns:  nothing.
 */
void stree_set_max_alloc(int size)
{
  if (size <= 0)
    max_alloc = 0;
  else
    max_alloc = size;
}






/*
 *
 *
 * Internal procedures to use when building and manipulating trees.
 *
 *
 *
 */


/*
 * int_stree_init_alphamaps
 *
 * Initialize the maps that map the raw stings into the 0..alphasize-1
 * strings.  There are three maps, a DNA/RNA map, a PROTEIN map and
 * and ASCII map.
 *
 * Parameters:  none
 * Returns:  nothing
 */
static void int_stree_init_alphamaps(void)
{
  int i;

  for (i=0; i < 128; i++)
    dnamap[i] = 4;
  dnamap['a'] = dnamap['A'] = 0;
  dnamap['c'] = dnamap['C'] = 1;
  dnamap['g'] = dnamap['G'] = 2;
  dnamap['t'] = dnamap['T'] = 3;
  dnamap['u'] = dnamap['U'] = 3;

  dnamap['r'] = dnamap['R'] = 0;
  dnamap['y'] = dnamap['Y'] = 1;
  dnamap['w'] = dnamap['W'] = 0;
  dnamap['s'] = dnamap['S'] = 1;
  dnamap['m'] = dnamap['M'] = 0;
  dnamap['k'] = dnamap['K'] = 2;
  dnamap['h'] = dnamap['H'] = 0;
  dnamap['b'] = dnamap['B'] = 1;
  dnamap['v'] = dnamap['V'] = 0;
  dnamap['d'] = dnamap['D'] = 0;
  dnamap['n'] = dnamap['N'] = 0;

  for (i=0; i < 128; i++)
    proteinmap[i] = 20;
  proteinmap['a'] = proteinmap['A'] = 0;
  proteinmap['c'] = proteinmap['C'] = 1;
  proteinmap['d'] = proteinmap['D'] = 2;
  proteinmap['e'] = proteinmap['E'] = 3;
  proteinmap['f'] = proteinmap['F'] = 4;
  proteinmap['g'] = proteinmap['G'] = 5;
  proteinmap['h'] = proteinmap['H'] = 6;
  proteinmap['i'] = proteinmap['I'] = 7;
  proteinmap['k'] = proteinmap['K'] = 8;
  proteinmap['l'] = proteinmap['L'] = 9;
  proteinmap['m'] = proteinmap['M'] = 10;
  proteinmap['n'] = proteinmap['N'] = 11;
  proteinmap['p'] = proteinmap['P'] = 12;
  proteinmap['q'] = proteinmap['Q'] = 13;
  proteinmap['r'] = proteinmap['R'] = 14;
  proteinmap['s'] = proteinmap['S'] = 15;
  proteinmap['t'] = proteinmap['T'] = 16;
  proteinmap['v'] = proteinmap['V'] = 17;
  proteinmap['w'] = proteinmap['W'] = 18;
  proteinmap['y'] = proteinmap['Y'] = 19;

  proteinmap['b'] = proteinmap['B'] = 2;
  proteinmap['z'] = proteinmap['Z'] = 3;
  proteinmap['x'] = proteinmap['X'] = 0;

  for (i=0; i < 128; i++)
    selfmap[i] = (char) i;
}


/*
 * int_stree_insert_string
 *
 * Insert a string into the list of strings maintained in the
 * SUFFIX_TREE structure, in preparation for adding the suffixes of
 * the string to the tree.
 *
 * Parameters:  tree  -  A suffix tree
 *              S     -  The sequence
 *              M     -  The sequence length
 *
 * Returns:  The internal index into the SUFFIX_TREE structure's
 *           strings/rawstrings/lengths/ids arrays, or -1 on an error.
 */
int int_stree_insert_string(SUFFIX_TREE tree, char *S, int M, int strid)
{

  int i, slot, next;
  char *buffer;

  if (tree->nextslot == MAXNUMSTR)
    return -1;

  for (i=0; i < M; i++)
	  //cout << "stree_mapch " << stree_mapch(tree,S[i]) << endl;
  	  //cout << "alpha " << tree->alpha_size << endl;
    if (stree_mapch(tree, S[i]) >= tree->alpha_size)
      return -1;

  if (tree->copyflag) {
    if ((buffer = (char*)malloc(M + 1)) == NULL)
      return -1;

    strncpy(buffer, S, M+1);
    S = buffer;
  }

  slot = tree->nextslot;
  tree->strings[slot] = S;
  tree->lengths[slot] = M;
  tree->ids[slot] = strid;

  for (next=slot+1; next < MAXNUMSTR; next++)
    if (tree->strings[next] == NULL)
      break;
  tree->nextslot = next;

  return slot;
}


/*
 * int_stree_delete_string
 *
 * Remove a string from the list of strings maintained in the
 * SUFFIX_TREE structure, as the last step in removing the
 * suffixes of that string from the tree.
 *
 * Parameters:  tree  -  A suffix tree
 *              id    -  The internal id for the string
 *
 * Returns: nothing
 */
void int_stree_delete_string(SUFFIX_TREE tree, int id)
{
  if (tree->strings[id] == NULL)
    return;

  if (tree->copyflag)
    free(tree->strings[id]);

  tree->strings[id] = NULL;
  if (id < tree->nextslot)
    tree->nextslot = id;
}


/*
 * int_stree_convert_leafnode
 *
 * Convert a LEAF structure into a NODE structure and replace the
 * NODE for the LEAF in the suffix tree..
 *
 * Parameters:  tree  -  a suffix tree
 *              node  -  a leaf of the tree
 *
 * Returns:  The NODE structure corresponding to the leaf, or NULL.
 */
STREE_NODE int_stree_convert_leafnode(SUFFIX_TREE tree, STREE_NODE node)
{
  STREE_NODE newnode;
  STREE_LEAF leaf;
  STREE_INTLEAF intleaf;

  if (!int_stree_isaleaf(tree, node))
    return node;

  leaf = (STREE_LEAF) node;

  newnode = int_stree_new_node(tree, stree_get_edgestr(tree, node),
                               stree_get_edgelen(tree, node));
  if (newnode == NULL)
    return NULL;

  intleaf = int_stree_new_intleaf(tree, leaf->strid,
                                  int_stree_get_leafpos(tree, leaf));
  if (intleaf == NULL) {
    int_stree_free_node(tree, newnode);
    return NULL;
  }

  newnode->id = leaf->id;
  newnode->isaleaf = 0;
  newnode->ch = 1;
  newnode->children = (STREE_NODE) intleaf;

  int_stree_reconnect(tree, stree_get_parent(tree, node), node, newnode);
  int_stree_free_leaf(tree, leaf);

  return newnode;
}


/*
 * int_stree_get_suffix_link
 *
 * Traverses the suffix link from a node, and returns the node at the
 * end of the suffix link.
 *
 * Parameters:  tree  -  a suffix tree
 *              node  -  a tree node
 *
 * Return:  The node at the end of the suffix line.
 */
STREE_NODE int_stree_get_suffix_link(SUFFIX_TREE tree, STREE_NODE node)
{
  int len, edgelen;
  char *edgestr;
  STREE_NODE parent;

  if (node == stree_get_root(tree))
    return NULL;
  else if (!int_stree_isaleaf(tree, node))
    return node->suffix_link;

  edgestr = stree_get_edgestr(tree, node);
  edgelen = stree_get_edgelen(tree, node);
  parent = stree_get_parent(tree, node);

  /*
   * Do the skip/count trip to skip down to the proper node.
   */
  if (parent != stree_get_root(tree))
    parent = parent->suffix_link;
  else {
    edgestr++;
    edgelen--;
  }

  node = parent;
  while (edgelen > 0) {
    node = stree_find_child(tree, node, *edgestr);
    assert(node != NULL);

    len = stree_get_edgelen(tree, node);
    edgestr += len;
    edgelen -= len;
  }

  return node;
}


/*
 * int_stree_get_parent
 *
 * Gets the parent of the current node.
 *
 * Parameters:  tree   -  A suffix tree
 *              node   -  A tree node
 *
 * Returns:  The parent node, or NULL if node is the root.
 */
STREE_NODE int_stree_get_parent(SUFFIX_TREE tree, STREE_NODE node)
{
  while (!node->nextisparent)
    node = node->next;

  return node->next;
}


/*
 * int_stree_get_leafpos
 *
 * Return the position of the suffix that ends at that leaf.  This value
 * must be reconstructed by walking back up the tree to the root.
 *
 * Parameters:  tree   -  A suffix tree
 *              node   -  A tree node
 *
 * Returns:  The edge label length.
 */
int int_stree_get_leafpos(SUFFIX_TREE tree, STREE_LEAF leaf)
{
  int pos;
  STREE_NODE node, root;

  pos = leaf->pos;

  node = (STREE_NODE) leaf;
  root = stree_get_root(tree);
  while (1) {
    while (!node->nextisparent)
      node = node->next;

    if (node == root)
      break;

    node = node->next;
    pos -= node->edgelen;
  }

  return pos;
}


/*
 * int_stree_connect
 *
 * Connect a node as the child of another node.
 *
 * Parameters:  tree   -  A suffix tree
 *              node   -  The node to get the new child.
 *              child  -  The child being added.
 *
 * Returns:  The parent after the child has been connected (if the
 *           parent was originally a leaf, this may mean replacing
 *           the leaf with a node).
 */
STREE_NODE int_stree_connect(SUFFIX_TREE tree, STREE_NODE parent,
                             STREE_NODE newchild)
{
  char ch;
  STREE_NODE node, back;

  if (int_stree_isaleaf(tree, parent) &&
      (parent = int_stree_convert_leafnode(tree, parent)) == NULL)
    return NULL;

  ch = stree_get_mapch(tree, newchild);
  node = stree_get_children(tree, parent);
  back = NULL;
  while (node != NULL && stree_get_mapch(tree, node) < ch) {
    back = node;
    node = stree_get_next(tree, node);
  }

  if (node != NULL) {
    if (stree_get_mapch(tree, node) == ch)
      return NULL;

    newchild->next = node;
    newchild->nextisparent = 0;
  }
  else {
    newchild->next = parent;
    newchild->nextisparent = 1;
  }

  if (back == NULL) {
    if (int_stree_has_intleaves(tree, parent))
      ((STREE_INTLEAF) (parent->children))->nextchild = newchild;
    else
      parent->children = newchild;
  }
  else {
    back->next = newchild;
    back->nextisparent = 0;
  }

  tree->idents_dirty = 1;
  return parent;
}


/*
 * int_stree_reconnect
 *
 * Replaces one node with another in the suffix tree, reconnecting
 * the link from the parent to the new node.
 *
 * Parameters:  tree      -  A suffix tree
 *              parent    -  The parent of the node being replaced
 *              oldchild  -  The child being replaced
 *              newchild  -  The new child
 *
 * Returns:  nothing
 */
int int_stree_reconnect(SUFFIX_TREE tree, STREE_NODE parent,
                        STREE_NODE oldchild, STREE_NODE newchild)
{
  STREE_NODE node, back;

  node = stree_get_children(tree, parent);
  back = NULL;
  while (node != NULL && node != oldchild) {
    back = node;
    node = stree_get_next(tree, node);
  }

  if (node == NULL)
    return 0;

  newchild->next = oldchild->next;
  newchild->nextisparent = oldchild->nextisparent;

  if (back == NULL) {
    if (int_stree_has_intleaves(tree, parent))
      ((STREE_INTLEAF) (parent->children))->nextchild = newchild;
    else
      parent->children = newchild;
  }
  else
    back->next = newchild;

  tree->idents_dirty = 1;
  return 1;
}


/*
 * int_stree_disc_from_parent
 *
 * Disconnect a node from its parent in the tree.
 * NOTE:  This procedure only does the link manipulation part of the
 *        disconnection process.  int_stree_disconnect is the real
 *        disconnection function.
 *
 * Parameters:  tree    -  A suffix tree
 *              parent  -  The parent node
 *              child   -  The child to be disconnected
 *
 * Return:  nothing.
 */
void int_stree_disc_from_parent(SUFFIX_TREE tree, STREE_NODE parent,
                                STREE_NODE child)
{
  STREE_NODE node, back;

  node = stree_get_children(tree, parent);
  back = NULL;
  while (node != NULL && node != child) {
    back = node;
    node = stree_get_next(tree, node);
  }

  if (node == NULL)
    return;

  node = stree_get_next(tree, node);
  if (back == NULL) {
    if (int_stree_has_intleaves(tree, parent))
      ((STREE_INTLEAF) (parent->children))->nextchild = node;
    else
      parent->children = node;
  }
  else {
    if (node == NULL) {
      back->next = parent;
      back->nextisparent = 1;
    }
    else
      back->next = node;
  }
}


/*
 * int_stree_disconnect
 *
 * Disconnects a node from its parent, and compacts the tree if that
 * parent is no longer needed.
 *
 * Parameters:  tree  -  a suffix tree
 *              node  -  a tree node
 *
 * Return:  The node at the end of the suffix line.
 */
void int_stree_disconnect(SUFFIX_TREE tree, STREE_NODE node)
{
  int num;
  STREE_NODE parent;

  if (node == stree_get_root(tree))
    return;

  parent = stree_get_parent(tree, node);
  int_stree_disc_from_parent(tree, parent, node);

  if (!int_stree_has_intleaves(tree, parent) &&
      parent != stree_get_root(tree) &&
      (num = stree_get_num_children(tree, parent)) < 2) {
    if (num == 0) {
      int_stree_disconnect(tree, parent);
      int_stree_delete_subtree(tree, parent);
    }
    else if (num == 1)
      int_stree_edge_merge(tree, parent);
  }

  tree->idents_dirty = 1;
}


/*
 * int_stree_edge_split
 *
 * Splits an edge of the suffix tree, and adds a new node between two
 * existing nodes at that split point.
 *
 * Parameters:  tree  -  a suffix tree
 *              node  -  The tree node just below the split.
 *              len   -  How far down node's edge label the split is.
 *
 * Return:  The new node added at the split.
 */
STREE_NODE int_stree_edge_split(SUFFIX_TREE tree, STREE_NODE node, int len)
{
  char *edgestr;
  STREE_NODE newnode, parent;
  STREE_LEAF leaf;

  if (tree->num_nodes == MAXNUMNODES || len == 0 ||
      stree_get_edgelen(tree, node) <= len)
    return NULL;

  edgestr = stree_get_edgestr(tree, node);

  newnode = int_stree_new_node(tree, edgestr, len);
  if (newnode == NULL)
    return NULL;

  parent = stree_get_parent(tree, node);
  int_stree_reconnect(tree, parent, node, newnode);

  if (int_stree_isaleaf(tree, node)) {
    leaf = (STREE_LEAF) node;
    leaf->pos += len;
    leaf->ch = stree_mapch(tree, edgestr[len]);
  }
  else {
    node->edgestr += len;
    node->edgelen -= len;
  }

  if (int_stree_connect(tree, newnode, node) == NULL) {
    if (int_stree_isaleaf(tree, node)) {
      leaf = (STREE_LEAF) node;
      leaf->pos -= len;
      leaf->ch = stree_mapch(tree, *edgestr);
    }
    else {
      node->edgestr -= len;
      node->edgelen += len;
    }
    int_stree_reconnect(tree, parent, newnode, node);
    int_stree_free_node(tree, newnode);
    return NULL;
  }

  tree->num_nodes++;
  tree->idents_dirty = 1;

  return newnode;
}


/*
 * int_stree_edge_merge
 *
 * When a node has no "leaves" and only one child, this function will
 * remove that node and merge the edges from parent to node and node
 * to child into a single edge from parent to child.
 *
 * Parameters:  tree  -  A suffix tree
 *              node  -  The tree node to be removed
 *
 * Return:  nothing.
 */
void int_stree_edge_merge(SUFFIX_TREE tree, STREE_NODE node)
{
  int len;
  STREE_NODE parent, child;
  STREE_LEAF leaf;

  if (node == stree_get_root(tree) || int_stree_isaleaf(tree, node) ||
      int_stree_has_intleaves(tree, node))
    return;

  parent = stree_get_parent(tree, node);
  child = stree_get_children(tree, node);
  if (stree_get_next(tree, child) != NULL)
    return;

  len = stree_get_edgelen(tree, node);
  if (int_stree_isaleaf(tree, child)) {
    leaf = (STREE_LEAF) child;
    leaf->pos -= len;
    leaf->ch = stree_get_mapch(tree, node);
  }
  else {
    child->edgestr -= len;
    child->edgelen += len;
  }

  int_stree_reconnect(tree, parent, node, child);
  tree->num_nodes--;
  tree->idents_dirty = 1;

  int_stree_free_node(tree, node);
}


/*
 * int_stree_add_intleaf
 *
 * Connects an intleaf to a node.
 *
 * Parameters:  tree  -  A suffix tree
 *              node  -  A tree node
 *              id    -  The internal identifier of the string
 *              pos   -  The position of the suffix in the string
 *
 * Returns:  Non-zero on success, zero on error.
 */
int int_stree_add_intleaf(SUFFIX_TREE tree, STREE_NODE node,
                          int strid, int pos)
{
  STREE_INTLEAF intleaf, next;

  if (int_stree_isaleaf(tree, node) ||
      (intleaf = int_stree_new_intleaf(tree, strid, pos)) == NULL)
    return 0;

  if (int_stree_has_intleaves(tree, node)) {
    next = (STREE_INTLEAF) node->children;
    intleaf->next = next;
    intleaf->nextchild = next->nextchild;
  }
  else {
    intleaf->nextchild = node->children;
  }

  node->children = (STREE_NODE) intleaf;
  node->ch = 1;
  return 1;
}


/*
 * int_stree_remove_intleaf
 *
 * Removes an intleaf from a node.
 *
 * Parameters:  tree  -  A suffix tree
 *              node  -  A tree node
 *              strid -  The internal identifier of the string
 *              pos   -  The position of the suffix in the string
 *
 * Returns:  Non-zero on success, zero on error.
 */
int int_stree_remove_intleaf(SUFFIX_TREE tree, STREE_NODE node,
                             int strid, int pos)
{
  STREE_INTLEAF intleaf, back;

  if (int_stree_isaleaf(tree, node) || !int_stree_has_intleaves(tree, node))
    return 0;

  back = NULL;
  intleaf = int_stree_get_intleaves(tree, node);
  for (back=NULL; intleaf != NULL; back=intleaf,intleaf=intleaf->next)
    if (intleaf->strid == strid && intleaf->pos == pos)
      break;

  if (intleaf == NULL)
    return 0;

  if (back != NULL)
    back->next = intleaf->next;
  else if (intleaf->next != NULL) {
    intleaf->next->nextchild = intleaf->nextchild;
    node->children = (STREE_NODE) intleaf->next;
  }
  else {
    node->children = intleaf->nextchild;
    node->ch = 0;
  }

  int_stree_free_intleaf(tree, intleaf);
  return 1;
}


/*
 * int_stree_delete_subtree
 *
 * Free up all of the memory associated with the subtree rooted at node.
 *
 * Parameters:  tree  -  a suffix tree
 *              node  -  a tree node
 *
 * Return:  nothing.
 */
void int_stree_delete_subtree(SUFFIX_TREE tree, STREE_NODE node)
{
  STREE_NODE child, next;
  STREE_INTLEAF intleaf, intnext;

  if (int_stree_isaleaf(tree, node))
    int_stree_free_leaf(tree, (STREE_LEAF) node);
  else {
    child = stree_get_children(tree, node);
    while (child != NULL) {
      next = stree_get_next(tree, child);
      int_stree_delete_subtree(tree, child);
      child = next;
    }

    if (int_stree_has_intleaves(tree, node)) {
      intleaf = (STREE_INTLEAF) node->children;
      while (intleaf != NULL) {
        intnext = intleaf->next;
        int_stree_free_intleaf(tree, intleaf);
        intleaf = intnext;
      }
    }

    int_stree_free_node(tree, node);
  }

  tree->idents_dirty = 1;
}


/*
 * int_stree_remove_to_position
 *
 * Remove the suffixes of a string from the suffix tree, and compact the
 * tree as necessary.
 *
 * NOTE:  This can be used either to remove a successfully added string,
 *        or a string which was only partially added to the tree (because
 *        an error stopped the add operation).  But it should only be used
 *        to completely remove a string.
 *
 * Parameters:  tree        -  A suffix tree
 *              id          -  The internally used id of the string to remove
 *              num_remove  -  How many positions to remove.
 *
 * Return:  nothing.
 */
void int_stree_remove_to_position(SUFFIX_TREE tree, int id, int num_remove)
{
  int M, walklen, rempos, pos, num, status;
  char *S;
  STREE_NODE node, next;
  STREE_LEAF leaf;

  if (num_remove == 0)
    return;

  S = int_stree_get_string(tree, id);
  M = int_stree_get_length(tree, id);

  walklen = int_stree_walk_to_leaf(tree, stree_get_root(tree), 0, S, M,
                                   &node, &pos);
  assert(walklen == M || int_stree_isaleaf(tree, node));

  next = NULL;
  rempos = 0;
  while (rempos < num_remove) {
    if (rempos < num_remove - 1) {
      next = stree_get_suffix_link(tree, node);
      assert(next != NULL);
    }

    if (int_stree_isaleaf(tree, node)) {
      leaf = (STREE_LEAF) node;
      assert(leaf->strid == id &&
             int_stree_get_leafpos(tree, leaf) == rempos);

      int_stree_disconnect(tree, node);
      int_stree_free_leaf(tree, leaf);
    }
    else {
      status = int_stree_remove_intleaf(tree, node, id, rempos);
      assert(status != 0);

      if (!int_stree_has_intleaves(tree, node) &&
          node != stree_get_root(tree) &&
          (num = stree_get_num_children(tree, node)) < 2) {
        if (num == 0) {
          int_stree_disconnect(tree, node);
          int_stree_delete_subtree(tree, node);
        }
        else if (num == 1)
          int_stree_edge_merge(tree, node);
      }
    }

    if (rempos < num_remove - 1)
      node = next;
    rempos++;
  }
}


/*
 * int_stree_walk_to_leaf
 *
 * Traverses the suffix link from a node, and returns the node at the
 * end of the suffix link.
 *
 * Parameters:  tree     -  A suffix tree
 *              node     -  The starting node of the walk
 *              pos      -  The starting point on the node's edge.
 *              T        -  The string to match
 *              N        -  The matching string's length
 *              node_out - Where the walk ends
 *              pos_out  - Where on the node's edge does the walk end.
 *
 * Return:  The number of characters matched during the walk.
 */
int int_stree_walk_to_leaf(SUFFIX_TREE tree, STREE_NODE node, int pos,
                           char *T, int N, STREE_NODE *node_out, int *pos_out)
{
  int len, edgelen;
  char *edgestr;
  STREE_NODE child;

  if (int_stree_isaleaf(tree, node)) {
    *node_out = node;
    *pos_out = pos;
    return 0;
  }

  edgestr = stree_get_edgestr(tree, node);
  edgelen = stree_get_edgelen(tree, node);
  len = 0;
  while (1) {
    while (len < N && pos < edgelen &&
           stree_mapch(tree, T[len]) == stree_mapch(tree, edgestr[pos])) {
      pos++;
      len++;
    }

    if (len == N || pos < edgelen ||
        (child = stree_find_child(tree, node, T[len])) == NULL)
      break;

    if (int_stree_isaleaf(tree, child)) {
      *node_out = child;
      *pos_out = 0;
      return len;
    }

    node = child;
    edgestr = stree_get_edgestr(tree, node);
    edgelen = stree_get_edgelen(tree, node);
    pos = 1;
    len++;
  }

  *node_out = node;
  *pos_out = pos;
  return len;
}


/*
 * int_stree_set_idents
 *
 * Uses the non-recursive traversal to set the identifiers for the current
 * nodes of the suffix tree.  The nodes are numbered in a depth-first
 * manner, beginning from the root and taking the nodes in the order they
 * appear in the children lists.
 *
 * Parameters:  tree  -  A suffix tree
 *
 * Return:  nothing.
 */
void int_stree_set_idents(SUFFIX_TREE tree)
{
  int id;
  STREE_NODE node, next;

  if (!tree->idents_dirty)
    return;

  tree->idents_dirty = 0;

  /*
   * Use a non-recursive traversal.  See stree_traverse_subtree for
   * details.
   */
  id = 0;
  node = stree_get_root(tree);
  while (1) {
    node->id = id++;

    next = stree_get_children(tree, node);
    if (next != NULL) {
      node = next;
      continue;
    }

    while (1) {
      if (node == stree_get_root(tree))
        return;
      if ((next = stree_get_next(tree, node)) != NULL)
        break;

      node = stree_get_parent(tree, node);
    }

    node = next;
  }
}


/*
 *
 *
 * The memory allocation functions.
 *
 *
 * A single page allocator is use to allocate pages of nodes, leaves and
 * intleaves.  The memory is allocated in pages of 27648 bytes (this is
 * a number must be a multiple of 12, 16 and 24 bytes).
 */

#define PAGESIZE 27648
#define NODES_PER_PAGE (PAGESIZE / sizeof(SNODE_STRUCT))
#define LEAFS_PER_PAGE (PAGESIZE / sizeof(SLEAF_STRUCT))
#define INTLEAFS_PER_PAGE (PAGESIZE / sizeof(SINTLEAF_STRUCT))

typedef enum { LEAF, INTLEAF, NODE } PAGETYPE;

typedef struct pagenode {
  void *start, *end;
  int count, size;
  struct pagenode *next;
} PAGE;

static PAGE *pagelist = NULL;


/*
 * getpage
 *
 * This returns a new page to the structure specific allocation function.
 *
 * Parameters:  type  -  The type of page allocated (used to determine
 *                       the number of structures per page).
 *
 * Returns: An allocated page, or NULL if too much memory has been
 *          allocated or malloc fails.
 */
static void *getpage(PAGETYPE type)
{
  void *buffer;
  PAGE *newpage;

  if (max_alloc > 0 && cur_alloc + PAGESIZE + sizeof(PAGE) > max_alloc)
    return NULL;

  if ((newpage = (PAGE*)malloc(sizeof(PAGE))) == NULL ||
      (buffer = (PAGE*)malloc(PAGESIZE)) == NULL) {
    if (newpage != NULL)
      free(newpage);
    return NULL;
  }

  memset(buffer, 0, PAGESIZE);

  newpage->start = buffer;
  newpage->end = buffer + PAGESIZE;
  newpage->count = 0;
  newpage->size = (type == NODE ? NODES_PER_PAGE
                                : (type == LEAF ? LEAFS_PER_PAGE
                                                : INTLEAFS_PER_PAGE));
  newpage->next = pagelist;
  pagelist = newpage;

  cur_alloc += PAGESIZE + sizeof(PAGE);

  return buffer;
}


/*
 * free_element
 *
 * "Frees" the memory for a structure by finding the page containing the
 * structure and incrementing a count of the number of structures in the
 * page that have been "freed".  When all of the structures in a page have
 * been "freed", the page itself is freed.
 *
 * Parameters:  address  -  The address of the structure to be "freed"
 *
 * Returns:  nothing
 */
static void free_element(void *address)
{
  PAGE *page, *back;

  back = NULL;
  for (page=pagelist; page != NULL; back=page,page=page->next)
    if (page->start <= address && address < page->end)
      break;

  assert(page != NULL);

  page->count++;
  if (page->count == page->size) {
    if (back == NULL)
      pagelist = pagelist->next;
    else
      back->next = page->next;

    free(page->start);
    free(page);

    cur_alloc -= PAGESIZE * sizeof(PAGE);
  }
}



/*
 * int_stree_new_intleaf
 *
 * Allocates memory for a new INTLEAF structure.
 *
 * Parameters:  strid  -  The id of the string containing the new suffix
 *                        to be added to the tree.
 *              pos    -  The position of the new suffix in the string.
 *
 * Returns:  The structure or NULL.
 */
STREE_INTLEAF int_stree_new_intleaf(SUFFIX_TREE tree, int strid, int pos)
{
  static int pagecount = 0;
  static STREE_INTLEAF page = NULL;
  STREE_INTLEAF intleaf;

  if (page == NULL && (page = (STREE_INTLEAF) getpage(INTLEAF)) == NULL)
    return NULL;

  intleaf = &page[pagecount];
  intleaf->strid = strid;
  intleaf->pos = pos;

  if (++pagecount == INTLEAFS_PER_PAGE) {
    page = NULL;
    pagecount = 0;
  }

  return intleaf;
}


/*
 * int_stree_new_leaf
 *
 * Allocates memory for a new LEAF structure.
 *
 * Parameters:  strid       -  The id of the string containing the new suffix
 *                             to be added to the tree.
 *              pos         -  The position of the new suffix in the string.
 *              edgestr     -  The edge label on the edge to the leaf.
 *              rawedgestr  -  The raw edge label on the edge to the leaf.
 *              edgelen     -  The edge label's length.
 *
 * Returns:  The structure or NULL.
 */
STREE_LEAF int_stree_new_leaf(SUFFIX_TREE tree, int strid, int edgepos)
{
  static int pagecount = 0;
  static STREE_LEAF page = NULL;
  STREE_LEAF leaf;

  if (page == NULL && (page = (STREE_LEAF) getpage(LEAF)) == NULL)
    return NULL;

  leaf = &page[pagecount];
  leaf->isaleaf = 1;
  leaf->ch = stree_mapch(tree, tree->strings[strid][edgepos]);
  leaf->strid = strid;
  leaf->pos = edgepos;

  if (++pagecount == LEAFS_PER_PAGE) {
    page = NULL;
    pagecount = 0;
  }

  return leaf;
}


/*
 * int_stree_new_node
 *
 * Allocates memory for a new NODE structure.
 *
 * Parameters:  edgestr     -  The edge label on the edge to the node.
 *              edgelen     -  The edge label's length.
 *
 * Returns:  The structure or NULL.
 */
STREE_NODE int_stree_new_node(SUFFIX_TREE tree, char *edgestr, int edgelen)
{
  static int pagecount = 0;
  static STREE_NODE page = NULL;
  STREE_NODE node;

  if (page == NULL && (page = (STREE_NODE) getpage(NODE)) == NULL)
    return NULL;

  node = &page[pagecount];
  node->edgestr = edgestr;
  node->edgelen = edgelen;

  if (++pagecount == NODES_PER_PAGE) {
    page = NULL;
    pagecount = 0;
  }

  return node;
}


/*
 * int_stree_free_{intleaf,leaf,node}
 *
 * Free the memory used for an INTLEAF, LEAF or NODE structure.  Also,
 * if the NODE structure uses a children array, free that space too.
 *
 * Parameters:  ileaf/leaf/node  -  The structure to free.
 *
 * Return:  nothing.
 */
void int_stree_free_intleaf(SUFFIX_TREE tree, STREE_INTLEAF intleaf)
{
  free_element(intleaf);
}

void int_stree_free_leaf(SUFFIX_TREE tree, STREE_LEAF leaf)
{
  free_element(leaf);
}

void int_stree_free_node(SUFFIX_TREE tree, STREE_NODE node)
{
  free_element(node);
}



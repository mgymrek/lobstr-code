#ifndef STREE_H
#define STREE_H

/*
 * NOTE: The values of NUMSTRBITS and STRLENBITS must add up to 32.
 *       At their current values of 13 and 19, the suffix tree data
 *       structure can handle 8192 strings, each of whose length can
 *       be up to 524,288 characters.
 *
 *       Also, the number of tree nodes is limited to be less than
 *       or equal to 8,388,608 nodes.  This is the maximum number of
 *       unique identifiers that can be put in the 23 bits allocated
 *       to the "id" field of the nodes.
 */
#define NUMSTRBITS 16
#define STRLENBITS 16

#define MAXNUMSTR (1 << NUMSTRBITS)
#define MAXSTRLEN (1 << STRLENBITS)
#define MAXNUMNODES (1 << 23)
#define MAXALPHA (1 << 7)

#define STREE_DNA -1
#define STREE_RNA -1
#define STREE_PROTEIN -2
#define STREE_ASCII 128

struct stree_node;

typedef struct stree_intleaf {
  unsigned int strid : NUMSTRBITS;
  unsigned int pos : STRLENBITS;

  struct stree_intleaf *next;
  struct stree_node *nextchild;
} SINTLEAF_STRUCT, *STREE_INTLEAF;


typedef struct stree_leaf {
  unsigned int id : 23;
  unsigned int isaleaf : 1;
  unsigned int nextisparent : 1;
  unsigned int ch : 7;

  struct stree_node *next;

  unsigned int strid : NUMSTRBITS;
  unsigned int pos : STRLENBITS;
} SLEAF_STRUCT, *STREE_LEAF;


typedef struct stree_node {
  unsigned int id : 23;
  unsigned int isaleaf : 1;
  unsigned int nextisparent : 1;
  unsigned int ch : 7;

  struct stree_node *next;
  struct stree_node *children;
  struct stree_node *suffix_link;

  char *edgestr;
  int edgelen;
} SNODE_STRUCT, *STREE_NODE;

typedef struct {
  STREE_NODE root;
  int num_nodes;

  char **strings;
  int *lengths, *ids;
  int nextslot, copyflag;

  int alpha_size, idents_dirty;
  char *alpha_map;
} STREE_STRUCT, *SUFFIX_TREE;


SUFFIX_TREE stree_new_tree(int alphabet, int copyflag);
void stree_delete_tree(SUFFIX_TREE tree);

int stree_add_string(SUFFIX_TREE tree, char *S, int M, int id);
int stree_remove_string(SUFFIX_TREE tree, int id);
int stree_match(SUFFIX_TREE tree, char *T, int N,
                STREE_NODE *node_out, int *pos_out);
int stree_walk(SUFFIX_TREE tree, STREE_NODE node, int pos, char *T, int N,
               STREE_NODE *node_out, int *pos_out);

STREE_NODE stree_find_child(SUFFIX_TREE tree, STREE_NODE node, char ch);
int stree_get_num_children(SUFFIX_TREE tree, STREE_NODE node);

#define stree_get_children(tree,node) \
  (int_stree_isaleaf(tree,node) || (node)->children == NULL	\
   ? NULL							\
   : (int_stree_has_intleaves(tree, node)			    \
      ? ((STREE_INTLEAF) ((node)->children))->nextchild		    \
      : (node)->children))
#define stree_get_next(tree,node)  ((node)->nextisparent ? NULL : (node)->next)

#define stree_get_root(tree) ((tree)->root)
#define stree_get_num_nodes(tree) ((tree)->num_nodes)
#define stree_get_parent(tree,node) \
    ((node)->nextisparent ? (node)->next : int_stree_get_parent(tree, node))
#define stree_get_suffix_link(tree,node)  \
    (!int_stree_isaleaf(tree,node) ? (node)->suffix_link \
     : int_stree_get_suffix_link(tree, node))
#define stree_get_edgestr(tree,node)		\
  (!int_stree_isaleaf(tree,node)		\
   ? (node)->edgestr					   \
   : ((tree)->strings[((STREE_LEAF) (node))->strid] +	   \
      ((STREE_LEAF) (node))->pos))
#define stree_get_edgelen(tree,node)		\
  (!int_stree_isaleaf(tree,node)		\
   ? (node)->edgelen					   \
   : ((tree)->lengths[((STREE_LEAF) (node))->strid] -	   \
      ((STREE_LEAF) (node))->pos))

// (dyoo) I added this macro to allow the user to recover the strid's
// they used to insert the strings.  Be careful not to confuse this with
// the internally used int_stree_get_strid() function.
#define stree_get_strid(tree, node) \
    (!int_stree_isaleaf(tree,node) \
        ? int_stree_get_strid(tree, (node)->id) \
        : int_stree_get_strid(tree, ((STREE_LEAF) (node))->strid))


/* Be careful: stree_getch() will segfault if node is the root node,
   so make sure not to call stree_getch() on root.  We may want to make
   an assertion to test this. */
#define stree_getch(tree,node) \
  (!int_stree_isaleaf(tree,node)		\
   ? ((node)->edgestr)[0]				    \
   : (((tree)->strings[((STREE_LEAF) (node))->strid] +	    \
       ((STREE_LEAF) (node))->pos))[0])

#define stree_mapch(tree,ch) ((tree)->alpha_map[(int) (ch)])
#define stree_get_mapch(tree,node) \
  (int_stree_isaleaf(tree,node)	   \
   ? (node)->ch						\
   : (tree)->alpha_map[(int) *((node)->edgestr)])


#define stree_get_ident(tree,node) \
  (!((tree)->idents_dirty) ? (node)->id					\
   : (int_stree_set_idents(tree), (node)->id))

int stree_get_labellen(SUFFIX_TREE tree, STREE_NODE node);
void stree_get_label(SUFFIX_TREE tree, STREE_NODE node, char *buffer,
                     int buflen, int endflag);

int stree_get_num_leaves(SUFFIX_TREE tree, STREE_NODE node);
int stree_get_leaf(SUFFIX_TREE tree, STREE_NODE node, int leafnum,
                   char **string_out, int *pos_out, int *id_out);

void stree_set_max_alloc(int size);


/*
 *
 * Internal procedures to use when building and manipulating trees.
 *
 */

int int_stree_insert_string(SUFFIX_TREE tree, char *S, int M, int strid);
void int_stree_delete_string(SUFFIX_TREE tree, int strid);

STREE_NODE int_stree_convert_leafnode(SUFFIX_TREE tree, STREE_NODE node);

#define int_stree_isaleaf(tree,node)  ((node)->isaleaf)
#define int_stree_has_intleaves(tree,node) \
              (int_stree_isaleaf(tree,node) ? 0 : (node)->ch)
#define int_stree_get_intleaves(tree,node) \
              (int_stree_has_intleaves(tree,node) \
                 ? ((STREE_INTLEAF) (node)->children) \
                 : NULL)

STREE_NODE int_stree_get_suffix_link(SUFFIX_TREE tree, STREE_NODE node);
STREE_NODE int_stree_get_parent(SUFFIX_TREE tree, STREE_NODE node);
int int_stree_get_leafpos(SUFFIX_TREE tree, STREE_LEAF leaf);

#define int_stree_get_string(tree,id)  ((tree)->strings[(id)])
#define int_stree_get_length(tree,id)  ((tree)->lengths[(id)])
#define int_stree_get_strid(tree,id)  ((tree)->ids[(id)])

STREE_NODE int_stree_connect(SUFFIX_TREE tree, STREE_NODE parent,
                             STREE_NODE newchild);
int int_stree_reconnect(SUFFIX_TREE tree, STREE_NODE parent,
                        STREE_NODE oldchild, STREE_NODE newchild);
void int_stree_disc_from_parent(SUFFIX_TREE tree, STREE_NODE parent,
                                STREE_NODE child);
void int_stree_disconnect(SUFFIX_TREE tree, STREE_NODE node);

STREE_NODE int_stree_edge_split(SUFFIX_TREE tree, STREE_NODE node, int len);
void int_stree_edge_merge(SUFFIX_TREE tree, STREE_NODE node);

int int_stree_add_intleaf(SUFFIX_TREE tree, STREE_NODE node,
                          int strid, int pos);
int int_stree_remove_intleaf(SUFFIX_TREE tree, STREE_NODE node,
                             int strid, int pos);

void int_stree_delete_subtree(SUFFIX_TREE tree, STREE_NODE node);
void int_stree_remove_to_position(SUFFIX_TREE tree, int strid, int num_remove);

int int_stree_walk_to_leaf(SUFFIX_TREE tree, STREE_NODE node, int pos,
                           char *T, int N, STREE_NODE *node_out, int *pos_out);

void int_stree_set_idents(SUFFIX_TREE tree);

STREE_INTLEAF int_stree_new_intleaf(SUFFIX_TREE tree, int strid, int pos);
STREE_LEAF int_stree_new_leaf(SUFFIX_TREE tree, int strid, int edgepos);
STREE_NODE int_stree_new_node(SUFFIX_TREE tree, char *edgestr, int edgelen);
void int_stree_free_intleaf(SUFFIX_TREE tree, STREE_INTLEAF ileaf);
void int_stree_free_leaf(SUFFIX_TREE tree, STREE_LEAF leaf);
void int_stree_free_node(SUFFIX_TREE tree, STREE_NODE node);

#endif

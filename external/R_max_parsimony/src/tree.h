#ifndef tree_h
#define tree_h

#include <stdbool.h>
// a struct that holds a hierarchical tree

// All dynamically allocated memory is from R; hence
// do no alloc or free anything in an h_tree
struct h_tree {
  int node_n;
  int leaf_n;
  // The edges are describe in two arrays
  // containing parents and children.
  int edge_n;
  int *edge_child;
  int *edge_parent;
  // leaf nodes have no entry in the parent array
  // The root node has no entry in the child array
  // as it is not the child of any other node.
  
  const unsigned char **leaf_states;
  int dim_n;
  unsigned int al_offset;
  unsigned int al_size;
  // each leaf node is described by 0 terminated char* array
  // of length var_n;
  // var_min and var_max give the max and mininum values that are allowed
  // var_min is used as an offset in a substitution matrix.
  // These can be created in R with intToUtf8() and back with
  bool is_good;
};



// A node in an hierarcical tree. This is created from an h_tree as described above
// the node_descriptions are taken from
// In a binary tree a node has up to 3 edges; Usually these are:
// parent x 1
// children x 2
// However, the root node in an unrooted tree (???) has three child connections
// and leaf nodes have only a single parent connection. 
struct ht_node {
  struct ht_node* edges[3];
  int node_i; 
  int edges_i[3];
  bool is_child[3]; 
  int *tree_lengths;   // signed so we can pass it back to R.
  int *inferred_state; // one for every dimension.. 
  int *state_delta;  // int because of R; signed is nice as well.
  // The child states from which the tree lengths were inferred
  // these will have the same dimensions as tree_lengths, but we
  // can use char as the substitution matrix and the leaf_states are encoded
  // as chars..
  // The child states are not that useful for two reasons; one is that we may
  // have many states with a minimal score; the other is that we may in fact
  // have three children (true for a root node). This makes it easir to redo
  // the analysis when going downstream.
  /* unsigned char *child_states_1; */
  /* unsigned char *child_states_2; */
  unsigned int edge_n;  // maximum of three, so wasteful to make an int
  bool is_leaf;  // 0 or 1
  bool length_determined;
  // the tree_lengths is a matrix that holds the associated sub-tree length
  // for each potential variable value. It has
  // var_n columns and al_size rows in row major order. 
};

// desc has to be 0-terminated, as that is the only way of getting
// a char array out of R.
int check_state(const unsigned char *desc, int al_offset, int al_size);
void check_tree( struct h_tree *tree );

struct h_tree make_tree(int *edge_child, int *edge_parent, int edge_n, int node_n, int leaf_n,
			int dim_n, const char **leaf_states, int al_offset, int al_size);

void free_tree(struct h_tree tree);

struct ht_node* make_nodes(struct h_tree *tree, int *sub_matrix, int *root_i);
void ht_nodes_free(struct ht_node *nodes, int l);
void sankoff_set_lengths(struct ht_node *node, int *sub_matrix, int al_size, int dim_n);
void sankoff_infer_states(struct ht_node *node, int *sub_matrix, int al_size, int dim_n);

#endif

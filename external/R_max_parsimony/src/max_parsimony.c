#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>  // for abs
#include "tree.h"

// A helper function to allow R to convert ints to simple ascii. This will 
// return a null terminated char vector for a given integer vector.

SEXP intToAscii(SEXP ints_r){
  if(TYPEOF(ints_r) != INTSXP || length(ints_r) < 1)
    error("vector must be integral and of positive length");
  size_t l = length(ints_r);
  int *ints = INTEGER(ints_r);
  SEXP ret_data = PROTECT(allocVector( STRSXP, 1 ));
  unsigned char *chars = malloc(sizeof(char) * (l + 1));
  chars[l] = 0;
  for(size_t i=0; i < l; ++i)
    chars[i] = (char)( abs(ints[i]) % 255 );
  // this will copy, and as such is wasteful, but I have not found a way to allocate
  // memory in R.
  SET_STRING_ELT( ret_data, 0, mkChar( (const char*)chars ) );
  free(chars);
  UNPROTECT(1);
  return(ret_data);
}

// tree_r        : a matrix with two columns: parent, child
//                 that describes connections within the tree
// tree_prop_r   : the number of nodes and leaf_nodes
// sub_matrix_r  : a substitution matrix.
// alphabet_r    : al_offset and al_size used for the substitution matrix
// leaf_states_r : a character vector, where each entry defines the
//                 state of a leaf_node.
//  Note that how missing values are handled depends on the substitution
//  tree and 
SEXP sankoff(SEXP tree_r, SEXP tree_props_r, SEXP sub_matrix_r, SEXP alphabet_r,
	     SEXP leaf_states_r)
{
  // apart from leaf_states_r, everything should be INTSXP. leaf_states_r should be
  // STRSXP
  if(TYPEOF(tree_r) != INTSXP || TYPEOF(tree_props_r) != INTSXP || TYPEOF(sub_matrix_r) != INTSXP
     || TYPEOF(alphabet_r) != INTSXP || TYPEOF(leaf_states_r) != STRSXP )
    error("incorrect data type, should be: int, int, int, int, character");
  SEXP tree_dims = getAttrib(tree_r, R_DimSymbol);
  if(length(tree_dims) != 2)
    error("the tree should be a two-dimensional matrix");
  int tree_nrow = INTEGER(tree_dims)[0];
  int tree_ncol = INTEGER(tree_dims)[1];
  if(tree_ncol != 2 || tree_nrow < 3)
    error("The tree matrix should have at least 3 rows and exactly two columns");
  int *tree_edges = INTEGER(tree_r);
  int *tree_parents = tree_edges;
  int *tree_children = tree_edges + tree_nrow;

  // tree_props_r should have length 2
  if(length(tree_props_r) != 2)
    error("tree_props should have a length of 2");
  int node_n = INTEGER(tree_props_r)[0];
  int leaf_n = INTEGER(tree_props_r)[1];
  if(leaf_n < 2 || node_n <= leaf_n)
    error("there should be more nodes than leaves");

  // sub_matrix should be a square matrix
  SEXP sub_matrix_dims_r = getAttrib(sub_matrix_r, R_DimSymbol);
  if(length(sub_matrix_dims_r) != 2)
    error("sub_matrix should be a matrix");
  int *sub_matrix_dims = INTEGER(sub_matrix_dims_r);
  if(sub_matrix_dims[0] != sub_matrix_dims[1])
    error("sub_matrix should be square");
  int *sub_matrix = INTEGER(sub_matrix_r);

  // al dims should be an integer vector of length 2
  if(length(alphabet_r) != 2)
    error("alphabet_r should have length 2");
  int al_offset = INTEGER(alphabet_r)[0];
  int al_size = INTEGER(alphabet_r)[1];

  if(al_size != sub_matrix_dims[0])
    error("the dimensions of the submatrix should be the same as the alphabet size");

  // leaf_states_r should have the same length as the number of leaf nodes specified
  if(length(leaf_states_r) != leaf_n)
    error("the length of leaf_states_r should be the same as leaf_n");

  int dim_n = length( STRING_ELT(leaf_states_r, 0));
  const char **leaf_states = malloc(sizeof(char*) * length(leaf_states_r));
  for(int i=0; i < length(leaf_states_r); ++i){
    SEXP s = STRING_ELT(leaf_states_r, i);
    if(length(s) != dim_n){
      free(leaf_states);
      error("all leaf_states must be the same length");
    }
    leaf_states[i] = CHAR( s );
  }
  // and we should now be able to make the tree and then infer it's ancestry...

  struct h_tree tree = make_tree( tree_children, tree_parents, tree_nrow, node_n, leaf_n,
				  dim_n, leaf_states, al_offset, al_size );
  int root_i = -1;
  struct ht_node* nodes = make_nodes( &tree, sub_matrix, &root_i);

  if(root_i >= 0){
    sankoff_set_lengths( nodes + root_i, sub_matrix, al_size, dim_n );
    Rprintf("Calling infer_states\n");
    sankoff_infer_states( nodes + root_i, sub_matrix, al_size, dim_n );
    Rprintf("infer_states returned\n");
  }
  
  // That should give us a tree with distances in all the roots. We haven't yet solved
  // the issue of what to do with the root that has three children but no parents.
  // But in any case it will be necessary to somehow validat the tree before we do anything
  // else.
  SEXP ret_data = PROTECT(allocVector( VECSXP, tree.node_n ));
  // each element of the vector will have three entries, edge, length and state
  const char *names[3] = {"edge", "length", "state"};
  SEXP names_r = PROTECT(allocVector(STRSXP, 3));
  for(int i=0; i < length(names_r); ++i)
    SET_STRING_ELT(names_r, i, mkChar(names[i]));
  
  for(int i=0; i < tree.node_n; ++i){
    // a matrix edges_i 3 columns * 2 (the indices of the edges, and is child)
    // a matrix of the tree_lengths al_size * dim_n (dim_n in columns)
    // we can add the child states later, but first let us see how this works.. 
    SET_VECTOR_ELT(ret_data, i, allocVector(VECSXP, 3));
    SEXP elt = VECTOR_ELT(ret_data, i);
    setAttrib( elt, R_NamesSymbol, names_r );
    // I think that allocMatrix should be safe with a zero size. But we'll find out
    // Edges connected to the node
    SET_VECTOR_ELT(elt, 0, allocMatrix( INTSXP, nodes[i].edge_n, 2 ));
    // The tree lenghts for the different states
    SET_VECTOR_ELT(elt, 1, allocMatrix( INTSXP, al_size, dim_n ));
    // The inferred states and the size of the shift from one state to the
    // other. (directional)
    SET_VECTOR_ELT(elt, 2, allocMatrix( INTSXP, dim_n, 2 ));
    int *node_edges = INTEGER(VECTOR_ELT(elt, 0));
    int *node_lengths = INTEGER(VECTOR_ELT(elt, 1));
    int *inferred_states = INTEGER(VECTOR_ELT(elt, 2));
    int *state_delta = inferred_states + dim_n;
    for(int j=0; j < nodes[i].edge_n; ++j){
      node_edges[j] = nodes[i].edges_i[j] + 1;  // possibly plus 1?
      node_edges[j + nodes[i].edge_n] = (int)nodes[i].is_child[j];
    }
    memcpy( node_lengths, nodes[i].tree_lengths, sizeof(int) * al_size * dim_n );
    memcpy( inferred_states, nodes[i].inferred_state, sizeof(int) * dim_n );
    memcpy( state_delta, nodes[i].state_delta, sizeof(int) * dim_n );
  }
  // But lets see if we can compile it first.
  free_tree(tree);
  ht_nodes_free( nodes, tree.node_n );

  // make a dummy data set..
  UNPROTECT(2);
  return(ret_data);
}

static const R_CallMethodDef callMethods[] = {
  {"intToAscii", (DL_FUNC)&intToAscii, 1},
  {"sankoff", (DL_FUNC)&sankoff, 5},
  {NULL, NULL, 0}
};

void R_init_max_parsimony(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

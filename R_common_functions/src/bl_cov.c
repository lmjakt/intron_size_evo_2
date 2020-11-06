#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <stdlib.h>

// In order to avoid unnecessary parsing of data.frame columns ask the user
// to provide the relevant columns as vectors of the appropriate type.
// That can easily be done in R and makes this coding easier.

int assert_type(SEXP arg, int type, int l){
  return( TYPEOF(arg) == type && length(arg) == l );
}

// extract number trailing last '_';
// if no _, returns -1
// no error checking for atoi; would be better to use
// strtol(), but a little more cumbersom.
int extract_i(const char *s){
  int pos = -1;
  for(int i=0; s[i]; ++i){
    if(s[i] == '_')
      pos = i + 1;
  }
  if(pos == -1)
    return(pos);
  return( atoi( s + pos ));
}  

// Note that this does not do as much error checking as the R_function, but it should
// be much, much, much faster
// q_lengths_r : the lenghs expected for the set of query sequences
//               this is necessary so that we can set up 0 coverage vectors
//               where blast did not find anything.
// bl_qlen_r   : the qlengths reported by blast. These should be the same
//               as those given by q_lengths_r, and mapped via the last 
//               int in the query id.
SEXP query_covs(SEXP q_lengths_r,
		  SEXP q_id_r, SEXP bl_q_len_r, SEXP q_start_r, SEXP q_end_r,
		  SEXP s_id_r, SEXP s_exclude_r)
{
  if(TYPEOF(q_lengths_r) != INTSXP || length(q_lengths_r) < 1)
    error("q_lengths_r should be an integer vector of positive length");
  int *q_lengths = INTEGER(q_lengths_r);
  int q_n = length(q_lengths_r);

  if(TYPEOF(q_id_r) != STRSXP || length(q_id_r) < 1)
    error("q_id_r should be a character vector with positive length");
  int al_n = length(q_id_r);
  if(!assert_type( bl_q_len_r, INTSXP, al_n) ||
     !assert_type( q_start_r, INTSXP, al_n) ||
     !assert_type( q_end_r, INTSXP, al_n) ||
     !assert_type( s_id_r, STRSXP, al_n) ||
     !assert_type( s_exclude_r, STRSXP, 1 ))
    error("invalid argument for one or more of: bl_q_len_r, q_start_r, q_end_r, s_id_r, s_exlude_r");
  
  int *bl_q_len = INTEGER(bl_q_len_r);
  int *q_start = INTEGER(q_start_r);
  int *q_end = INTEGER(q_end_r);

  int **cov = malloc(sizeof(int*) * q_n );

  SEXP ret_val = PROTECT( allocVector( VECSXP, q_n ));
  for(int i=0; i < q_n; ++i){
    if( q_lengths[i] == NA_INTEGER || q_lengths[i] < 1 )
      continue;
    SET_VECTOR_ELT( ret_val, i, allocVector(INTSXP, q_lengths[i]) );
    cov[i] = INTEGER( VECTOR_ELT(ret_val, i) );
    memset( (void*)cov[i], 0, sizeof(int) * q_lengths[i] );
  }

  // and then simply go through and check...
  const char *s_exclude = CHAR(STRING_ELT( s_exclude_r, 0));
  
  for(int i=0; i < al_n; ++i){
    const char *q_id = CHAR(STRING_ELT( q_id_r, i));
    const char *s_id = CHAR(STRING_ELT( s_id_r, i));
    if( strstr( s_id, s_exclude ))
      continue;
    int j = extract_i( q_id ) - 1;
    if(j < 0 || j >= q_n){
      warning("unexpected query identifier: %s -> %d\n", q_id, j);
      continue;
    }
    if(q_lengths[j] != bl_q_len[i]){
      warning("non-matching query lengths: %d != %d", q_lengths[j], bl_q_len);
      continue;
    }
    int start = q_start[i] - 1;
    int end = q_end[i];
    if(start < 0 || end > q_lengths[j] || end < start){
      warning("obtained unreasonable query coordinates: %d -> %d (%d)", start, end, q_lengths[j]);
      continue;
    }
    for(int k=start; k < end; ++k)
      cov[j][k]++;
  }
  
  // and that should be everything except remember to clear the allocated memory
  free( cov );
  UNPROTECT(1);
  return( ret_val );
}

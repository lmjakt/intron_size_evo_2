#ifndef _CIGAR_H
#define _CIGAR_H

#include <stdint.h>
#include <string.h>

/*
  A cigar string defined by an alignment of a short transcript to a long
  genomic locus. The alignment is defined by the following operations
  0: aligned (two nucleotides are aligned)
  1: insertion into the transcript (gap character in the genomic locus)
  2: deletion from the transcript (gap character in the transcript)
  3: intron (gap character in the transcript, but matches to I positions
     in the transcript sequence)

  Note that the structs and functions defined here are pretty similar to those
  defined in samtools.h. These functions though are probably somewhat more
  limited and I am redefining them because my needs are quite specific.
 */

// According to the standard it is not necessary to specify the values
// as the first value will be 0, with subsequence values increasing by
// one if nothing is written. However, by making the values explicit here
// it is clear that the op part requires two bits;
typedef enum {
  MATCH=0,
  INSERT=1,
  DELETE=2,
  INTRON=3
} cigar_operation;

// Store the cigar string in unsigned 32 bit integers (uint32_t)
// Store the cigar_operation in the least significant bits
static const unsigned int op_mask = 3;
static const unsigned int count_shift = 2;

typedef struct {
  uint32_t *cigar;
  size_t capacity;
  size_t size;
} cigar_string;

// this cannot check if the string exists or not
void init_cigar(cigar_string *cigar, size_t capacity){
  cigar->capacity = capacity;
  cigar->size = 0;
  if(capacity){
    cigar->cigar = malloc( sizeof(uint32_t) * capacity );
    memset((void*)cigar->cigar, 0, sizeof(uint32_t) * capacity );
  }else{
    capacity = 0;
  }
}

void grow_cigar(cigar_string *cigar){
  size_t new_cap = cigar->capacity > 0 ? cigar->capacity * 2 : 6;
  uint32_t *old_cigar = cigar->cigar;
  cigar->cigar = malloc( sizeof(uint32_t) * new_cap );
  memset( (void*)cigar->cigar, 0, sizeof(uint32_t) * new_cap );
  memcpy((void*)cigar->cigar, old_cigar, cigar->size * sizeof(uint32_t));
  cigar->capacity = new_cap;
  free( old_cigar );
}

void push_cigar(cigar_string *cigar, cigar_operation op){
  if(cigar->size && cigar->cigar[ cigar->size-1 ] & op_mask == op){
    // then increment the counter and return;
    // this could do with some temporary pointers for clarity
    cigar->cigar[ cigar->size-1 ] = 
      ((1 + cigar->cigar[ cigar->size-1 ] >> count_shift) << count_shift) | op;
    return;
  }
  if(cigar->size >= cigar->capacity)
    grow_cigar(cigar);
  ++cigar->size;
  cigar->cigar[ cigar->size-1 ] = (1 << count_shift) | op;
}

void clear_cigar(cigar_string *cigar){
  memset( (void*)cigar->cigar, 0, sizeof(uint32_t) * cigar->capacity );
  cigar->size = 0;
}

void free_cigar(cigar_string *cigar){
  free(cigar->cigar);
  cigar->cigar = 0;
}

// convenience function to extract values

uint32_t cigar_value(cigar_string *cigar, size_t i){
  if(i < cigar->size )
    return( cigar->cigar[i] >> count_shift);
  return(0);
}

cigar_operation cigar_op(cigar_string *cigar, size_t i){
  if( i < cigar->size )
    return( cigar->cigar[i] & op_mask );
  return(0);
}

void alignment_lengths(cigar_string *cigar, 
		       uint32_t *tr_length, uint32_t *loc_length, uint32_t *al_length){
  tr_length = 0;
  loc_length = 0;
  al_length = 0;
  for(size_t i=0; i < cigar->size; ++i){
    uint32_t op_count = cigar->cigar[i] >> count_shift;
    switch( cigar->cigar[i] & op_mask ){
    case MATCH:
      tr_length += op_count;
      loc_length += op_count;
      al_length += op_count;
      break;
    case INSERT:
      tr_length += op_count;
      al_length += op_count;
      break;
    case DELETE:
      loc_length += op_count;
      al_length += op_count;
      break;
    case INTRON:
      loc_length += op_count;
      al_length += op_count;
      tr_length += 1;  // special case for this purpos;
    default:
      // do nothing;
    }
  }
}

#endif

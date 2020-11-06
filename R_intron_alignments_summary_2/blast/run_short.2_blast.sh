#!/bin/bash

## for the control sequences obtained after the non-control sequences
##
## Correction due to errors returned for mclapply when aligning control.control
## sequences in ../../R_intron_alignments_6/

## Run blast with $1 as the query file againts the databases listed in
## blast_db_list (created using find_dbs.sh

## the basic command without db and query

db_names=blast_db_list
queries=../extracted_seqs/short.2_*.fasta

for query in $queries
do
    cat blast_db_list | xargs -P 1 -I{} ./run_single_blast_2.sh $query {}
done



#!/bin/bash

## for the control sequences obtained after the non-control sequences

## Run blast with $1 as the query file againts the databases listed in
## blast_db_list (created using find_dbs.sh

## the basic command without db and query

db_names=blast_db_list
queries=../extracted_seqs/*.ctl_*.fasta

for query in $queries
do
    cat blast_db_list | xargs -I{} ./run_single_blast.sh $query {}
done



#!/bin/bash

## This runs blast against sequences in long_2 files. Note that the first
## run_blast.sh would also do this, except for that fact that those files
## did not exist at the time of writing.

## Run blast with $1 as the query file againts the databases listed in
## blast_db_list (created using find_dbs.sh

## the basic command without db and query

db_names=blast_db_list
queries=../*long_2.fasta

for query in $queries
do
    cat blast_db_list | xargs -I{} ./run_single_blast.sh $query {}
done



#!/bin/bash

## as for run_blast.sh, but only for the med files as these were
## added later on.

db_names=blast_db_list
queries=../extracted_seqs/*med*.fasta

for query in $queries
do
    cat blast_db_list | xargs -I{} ./run_single_blast.sh $query {}
done



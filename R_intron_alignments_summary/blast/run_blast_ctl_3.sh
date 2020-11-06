#!/bin/bash

## as for run_blast.sh and run_blast_med.sh, 
## Redone as I forgot to export mammalian control sequences
## in the wrong file.
##

db_names=blast_db_list
queries=../extracted_seqs/*.ctl*_mammalia.fasta

for query in $queries
do
    cat blast_db_list | xargs -I{} ./run_single_blast.sh $query {}
done



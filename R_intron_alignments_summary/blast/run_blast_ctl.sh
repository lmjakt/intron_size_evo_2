#!/bin/bash

## as for run_blast.sh and run_blast_med.sh, 
## but redone for the control set of introns as I found that I had read
## in the wrong file.

db_names=blast_db_list
queries=`(ls ../extracted_seqs/*ctl*fasta) | grep -v "\.ctl"`

for query in $queries
do
    cat blast_db_list | xargs -I{} ./run_single_blast.sh $query {}
done



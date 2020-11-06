#!/bin/bash

## as for run_blast.sh and run_blast_med.sh, 
## but redone for the long.ctl and ctl.ctl sequences as it turns out I
## exported the wrong sequences. 
## in the wrong file.
##

db_names=blast_db_list
queries=../extracted_seqs/*.ctl*.fasta

for query in $queries
do
    cat blast_db_list | xargs -I{} ./run_single_blast.sh $query {}
done



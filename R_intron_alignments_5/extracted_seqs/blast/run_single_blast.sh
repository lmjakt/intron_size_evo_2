#!/bin/bash

query=$1
db=$2

of=`basename $query .fasta`_`basename $db`.bl

blastn -query $query -db $db -out $of -task blastn -evalue 1e-5 -use_index true \
    -max_target_seqs 50 -num_threads 12 \
    -outfmt '6 qseqid sseqid qlen qstart qend sstart send evalue bitscore score length pident nident qcovs qcovhsp'



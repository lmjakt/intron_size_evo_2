#!/bin/bash

script=$HOME/genomes/ncbi_taxonomy/perl/is_member_of_taxon.pl
names=$HOME/genomes/ncbi_taxonomy/names.dmp
nodes=$HOME/genomes/ncbi_taxonomy/nodes.dmp

cut -f 2 Chordata_class_sizes.tsv | $script $names $nodes Vertebrata "clade" > vertebrate_membership.tsv

cut -f 2 Actinopterygii_order_sizes.tsv | $script $names $nodes Teleostei "infraclass" > teleost_membership.tsv

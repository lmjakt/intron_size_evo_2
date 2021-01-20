#!/bin/bash

## Get a table of chordate class sizes:
script=$HOME/genomes/catalogue_of_life/perl_scripts/species_counts.pl
tax_file=$HOME/genomes/catalogue_of_life/2020_12_01/Taxon.tsv
prof_file=$HOME/genomes/catalogue_of_life/2020_12_01/SpeciesProfile.tsv

echo $tax_file
echo $prof_file

$script $tax_file $prof_file PHYLUM Chordata 1 > Chordata_class_sizes.tsv
$script $tax_file $prof_file CLASS Actinopterygii 1 > Actinopterygii_order_sizes.tsv


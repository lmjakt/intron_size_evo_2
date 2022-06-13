#!/usr/bin/perl -w

## read in files in
## protein_align_features/

@al_files = glob( "protein_align_features/*feature.txt" );
@tax_files = @al_files;
map( {$_ =~ s/feature\.txt$/feature_tax.txt/;} @tax_files );


print join("\n", @al_files), "\n";
print "\n", join("\n", @tax_files), "\n";

## uniprot_file;
$uniprot_file = $ENV{HOME}."/genomes/uniprot/idmapping_selected.tab.gz";

## taxon file;
$tax_names = $ENV{HOME}."/genomes/ncbi_taxonomy_2021_09/names.dmp";

## first lets read in the names so we can get an id to string mapping:
open(IN, $tax_names ) || die "unable to open $tax_names $!\n";
while(<IN>){
    chomp;
    $_ =~ s/\|$//;
    @tmp = split /\t\|\t/, $_;
    $tmp[3] =~ s/\s+$//;
    if($tmp[3] eq "scientific name"){
	$sci_name{$tmp[0]} = $tmp[1];
    }
}
close(IN);

## Read in the uniprot file first. This takes something like 30 minutes
## to do. Most of the time seems to be related to the decompression.
## Note that the accession ids here do not have version numbers so
## we need to remove those from the others

open(IN, "-|", "zcat $uniprot_file") || die "unable to zcat $uniprot_file $!\n";
while(<IN>){
    @tmp = split /\t/, $_;
    $uni_tax{$tmp[0]} = $tmp[12];
}

## then we go through each file containing protein align information:
## there is one header row to discard;
## and then we have the following columns:
## 0. row
## 1. gene
## 2. transcript
## 3. hit_name
## 4. evalue
## 5. perc_ident
## 6. db_name
##
## We can map only if the db_name matches to UniProt

for $i(0..$#al_files){
    open(IN, $al_files[$i]) || die "unable to open $al_files[$i] for reading $!\n";
    open(OUT, ">$tax_files[$i]") || die "unable to open $tax_files[$i] for writing $!\n";
    chomp($header = <IN>);
    print OUT $header, "\ttax.id\ttaxon\n";
    while(<IN>){
	chomp;
	print OUT $_;
	@tmp = split /\t/, $_;
	if( $tmp[6] !~ /Uniprot/i ){
	    print OUT "\tNA\tNA\n";
	    next;
	}
	$prot_id = $tmp[3];
	$prot_id =~ s/\.\d+$//;
	$tax = $uni_tax{$prot_id} // "NA";
	$taxon = $sci_name{$tax} // "NA";
	print OUT "\t$tax\t$taxon\n";
    }
    close(OUT);
    close(IN);
}
	




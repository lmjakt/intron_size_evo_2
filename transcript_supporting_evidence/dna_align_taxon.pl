#!/usr/bin/perl -w

## read in files in
## protein_align_features/

@al_files = glob( "dna_align_features/*feature.txt" );
@tax_files = @al_files;
map( {$_ =~ s/feature\.txt$/feature_tax.txt/;} @tax_files );


print join("\n", @al_files), "\n";
print "\n", join("\n", @tax_files), "\n";

## accession id to taxon
$accession_file = $ENV{HOME}."/genomes/ncbi_taxonomy_2021_09/nucl_gb.accession2taxid.gz";

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

## 
## Read in the accession to taxonomy id. Ignore version number

open(IN, "-|", "zcat $accession_file") || die "unable to zcat $accession_file $!\n";
while(<IN>){
    @tmp = split /\t/, $_;
    $acc_tax{$tmp[0]} = $tmp[2];
}

## then we go through each file containing DNA align information:
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
## We can map only if the db_name EMBL or genbank

for $i(0..$#al_files){
    open(IN, $al_files[$i]) || die "unable to open $al_files[$i] for reading $!\n";
    open(OUT, ">$tax_files[$i]") || die "unable to open $tax_files[$i] for writing $!\n";
    chomp($header = <IN>);
    print OUT $header, "\ttax.id\ttaxon\n";
    while(<IN>){
	chomp;
	print OUT $_;
	@tmp = split /\t/, $_;
	$hit_id = $tmp[3];
	$hit_id =~ s/\.\d+$//;
	$tax = $acc_tax{$hit_id} // "NA";
	$taxon = $sci_name{$tax} // "NA";
	print OUT "\t$tax\t$taxon\n";
    }
    close(OUT);
    close(IN);
}
	




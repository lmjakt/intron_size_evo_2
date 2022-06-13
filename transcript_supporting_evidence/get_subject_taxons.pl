#!/bin/perl -w

## Heavily rewritten by Martin.
## Last change: restrict names to scientific names only if available.

$bl_file = "../miseq_random_blast/random_miseq_seqs_vs_nt2.bl";
@acc_map = ("taxonomy/dead_nucl.accession2taxid", "taxonomy/nucl_gb.accession2taxid");
$names_file = "taxonomy/names.dmp";
$nodes_file = "taxonomy/nodes.dmp";

# $acc_map = "/DATA/MarineGenomics/ama/data/taxonomy/nucl_gb.accession2taxid";
# $names_file = "/DATA/MarineGenomics/ama/data/taxonomy/names.dmp";
# $nodes_file = "/DATA/MarineGenomics/ama/data/taxonomy/nodes.dmp";

## read in the taxonomy tree
open($fh2, "<", $nodes_file) or die "unable to open $nodes_file $!\n";
while(<$fh2>){
    chomp;
    $_ =~ s/\s+\|\s+/|/g;
    @tmp = split /\|/, $_;
    for $i(0..$#tmp){
        $tmp[$i] =~ s/^\s+//;
        $tmp[$i] =~ s/\s+$//;
	$tax_group{$tmp[0]} = $tmp[2];
	$parent{$tmp[0]} = $tmp[1];
	## if you want to be able to trace the tree
	## in the opposite direction then you can also
	## do
	$daughters{$tmp[1]}{$tmp[0]} = 1;
    }
}
close($fh2);
print STDERR "Finished reading in the tree structure\n";


## read in the names
open(IN, $names_file) || die "unable to open $names_file $!\n";
while(<IN>){
    @tmp = split /\s*\|\s*/, $_;
    if($tmp[3] eq 'scientific name' || !defined($names{$tmp[0]}) ){
	$names{$tmp[0]} = $tmp[1];
    }
}
$names{-1} = "NA";

print STDERR "Read in the names\n";

## then read in the acc_map..
## use the simpler format of the data point.. 
for $acc_map_f(@acc_map){
    open(IN, $acc_map_f) || die "unable to open $acc_map $!\n";
    $tmp = <IN>; ## get rid of the header
    while(<IN>){
	@tmp = split /\t/, $_;
	$acc_tax{$tmp[0]} = $tmp[2];
    }
    print STDERR "Read from : $acc_map_f\n";
}
print STDERR "Read in the accession to tax ids\n";


## The blast file is much, much smaller than the accession2id file so we go
## through the blast file first and make a hash of alignments with the blast
## data as an array;
open(IN, $bl_file) || die "unable to open $bl_file $!\n";
## only analyse the first alignment for each query
%observed_queries = ();
while(<IN>){
    @tmp = split /\t/, $_;
    $query_id = $tmp[0];
    next if(defined($observed_queries{$query_id}));
    $observed_queries{$query_id} = 1;
    if($tmp[1] =~ /\w+\|(\w+)\.?\d*\|/ || $tmp[1] =~ /(\w+)\.?\d*/){
	$id = $1;
    }else{
	print STDERR "unexpected identifier in $tmp[1]\n";
	next;
    }
    if(!defined($acc_tax{$id})){
    	print STDERR "No taxon for $id\n";
    	next;
    }
    $tax_id = $acc_tax{$id};
    @taxons = get_taxons($tax_id);
    @taxon_names = map( $names{$_}, @taxons );
    ## we simply print to STDOUT
    $tax_name = defined($names{$tax_id}) ? $names{$tax_id} : "NA";
    print STDOUT $id, "\t", $tax_id, "\t", $tax_name, "\t", ## $kingdom, "\t", $names{$kingdom}, "\t";
    join(",", @taxons), "\t", join(",", @taxon_names), "\t";
    print STDOUT $tmp[0], "\t", $tmp[12], "\t", $tmp[13], "\t", $tmp[14], "\t", $tmp[15], "\t", $tmp[16], "\n";
}

print STDERR "Read in the blast data\n";




## then if you want to retrace from daughter to the root node ## from C for example:
sub get_taxons {
    my $id = shift @_;
    my $node = $id;
    my @taxons = ($id);
    if(!defined($parent{$node})){
	print STDERR "No parent defined for $id\n";
	return(-1);
    }
    while( $tax_group{$node} ne "superkingdom" && $node != $parent{$node} ){
	$node = $parent{$node};
	push @taxons, $node;
    }
    return(@taxons);
}

## whilst doing something useful with the $node values (eg. sticking them ## in an array).
    
## to go from the root node and visit all other nodes in the tree you (probably) ## need to use recursion through the $daughters structure.
## something like:
    
sub traverse_tree {
    my $node = shift @_;
    ## do something useful with $node
    return if !defined($daughters{$node}) ;
    for my $daughter(keys %{$daughters{$node}} ){
	traverse_tree($daughter);
    }
}
    
    

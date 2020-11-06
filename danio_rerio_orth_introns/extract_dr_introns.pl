#!/usr/bin/perl -w

## extract the Danio rerio introns that have entries in the orthologue table
## dump these to a single fasta file with identifiers that will allow us to make some sense.

$fam_file = "../R_172_genomes/dr_intron_fam.txt";
$i_file = "../R_172_genomes/dr_intron_orthology_i.txt";
##$l_file = "../R_172_genomes/dr_intron_orthology_i.txt";
$tr_file = "../R_172_genomes/dr_intron_orthology_tr.txt";
$id_file = "../R_172_genomes/dr_intron_orthology_id.txt";

$column_name = 'danio rerio';

## read in the family file
open(IN, $fam_file) || die "unable to open $fam_file $!\n";
chomp($header = <IN>);
while(<IN>){
    chomp;
    @tmp = split /\t/, $_;
    push @fams, $tmp[1];
}

print STDERR "reading in orthology information\n";
@orth_i = extract_column($i_file, $column_name);
@orth_tr = extract_column($tr_file, $column_name);
@orth_id = extract_column($id_file, $column_name);


$fam_id = "";
%seq = ();

for $i(0..$#fams){
    if($fam_id ne $fams[$i]){
	$fam_id = $fams[$i];
	print STDERR "Reading sequence for $fam_id\n";
	%seq = read_fasta("../family_members/orthologue_transcripts/introns_".$fam_id)
    }
    $tr = $orth_tr[$i];
    $id = $orth_id[$i];
    $int_i = $orth_i[$i];
    die("no sequence defined for $tr\n") if !defined($seq{$tr});
    
    print ">", $id, "_", $tr, "_", $int_i, "_", $i, "\n", $seq{$tr}{seq}[$int_i-1], "\n";
}


sub extract_column {
    my($file, $column) = @_;
    open(IN, $file) || die "unable to open $i_file $!\n";
    chomp(my $header = <IN>);
    my @sp = split /\t/, $header;
    my $col_i = -1;
    for(my $i=0; $i < @sp; $i++){
	if($sp[$i] eq $column){
	    $col_i = $i;
	    last;
	}
    }
    die "unable to find $column\n" if $col_i == -1;
    my @data = ();
    while(<IN>){
	chomp;
	my @tmp = split /\t/, $_;
	push @data, $tmp[ $col_i ];
    }
    return(@data);
}

## a somewhat specialised read_fasta file.. 
sub read_fasta {
    my($file) = @_;
    open(IN, $file) || die "unable to open $file $!\n";
    my $gene_id = "";
    my $tr_id = "";
    my %seq = ();
    while(<IN>){
	chomp;
	if( $_ =~ /^>(\S+)/ ){
	    $gene_id = $1;
	    my @tmp = split /\t/, $_;
	    $tr_id = $tmp[1];
	    $seq{$tr_id}{gene} = $gene_id;
	    @{$seq{$tr_id}{seq}} = ();
	    next;
	}
	push @{$seq{$tr_id}{seq}}, $_;
    }
    return(%seq);
}
    
	    

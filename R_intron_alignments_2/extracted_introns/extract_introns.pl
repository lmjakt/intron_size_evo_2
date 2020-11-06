#!/usr/bin/perl -w

$f = $ARGV[0];
$of = $f;
$of =~ s/\.txt/.fasta/;

$prefix = "../../family_members/orthologue_transcripts/introns_";

open(OUT, ">", $of) || die "unable to open $of $!\n";

open(IN, $f) || die "unable to open $f $!\n";
while(<IN>){
    chomp;
    ($fam, $gene, $transcript, $intron, $orth_row) = split;
    $sfile = $prefix.$fam;
    open(SEQ, $sfile) || die "unable to open $sfile $!\n";
    while(<SEQ>){
	chomp;
	if($_ =~ /^>/){
	    @tmp = split /\t/, $_;
	    last if ($tmp[1] eq $transcript);
	}
    }
    $seq = "";
    for($i=0; $i < $intron; $i++){
	chomp($seq = <SEQ>);
    }
    ## And this should be the target intron. 
    print OUT ">$transcript", "_", $intron, "\t", $gene, "\t", 
    $transcript, "\t", $intron, "\t", $fam, "\t", $orth_row, "\n", $seq, "\n";
}


	
    

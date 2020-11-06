#!/usr/bin/perl -w

## open the specified fasta files and count the number of different residues in them

for $infile(@ARGV){
    open(IN, $infile) || die "unable to open $infile $!\n";
    print STDERR ".";
    while(<IN>){
	chomp;
	next if($_ =~ /^>/);
	for($i=0; $i < length($_); ++$i){
	    $counts{ substr($_, $i, 1) }++;
	}
    }
}
print STDERR "\n";

for $k( sort keys %counts ){
    print "$k\t$counts{$k}\n";
}

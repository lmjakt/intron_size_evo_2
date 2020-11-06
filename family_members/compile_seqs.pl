#!/usr/bin/perl -w

## call this with something like
## ./compile_seqs.pl vertebrate_family_members_1_15.txt genome_list

## compile RNA and protein sequences for orthologue sets from a set of different databases
## do most of the logic in perl; this means that we may end up loading rather more than we need.
## well lets see
$infile = shift @ARGV;
$seq_file_list = shift @ARGV;
use DBI;

## To connect to ensembl databases use the following pattern where we simply change the database name.
## $db = DBI->connect("DBI:mysql:database=ensembl_compara_97:host=ensembldb.ensembl.org:port=5306", "anonymous", "") || die $DBI::errstr;
## $db->disconnect();

## first get the names of the sequence files
open(IN, $seq_file_list) || die "unable to open $seq_file_list $!\n";
while(<IN>){
    chomp;
    push @seq_files, $_;
}

open(IN, $infile) || die "unable to open $infile $!\n";
chomp($line = <IN>);
@db_names = split /\t/, $line;
for $i(0..$#db_names){
    $sp_names[$i] = $db_names[$i];
    $sp_names[$i] =~ s/(^[^_]+_[^_]+)_.+/$1/;
##    print STDERR "db_map $sp_names[$i]  => $db_names[$i]\n";
    $db_map{$sp_names[$i]} = $db_names[$i];
}


while(<IN>){
    chomp;
    @tmp = split /\t/, $_, -1;
    if(@tmp != (@db_names + 1)){
	die "unexepted number of columns in \n$_\n";
    }
    $fam_id = $tmp[0];
    for $i(1..$#tmp){
	next if($tmp[$i] !~ /\d+/);
	@ids = split /,/, $tmp[$i];
	for $id(@ids){
	    $gene_ids{$db_names[$i-1]}{$id} = 1;
	    $fam{$fam_id}{$id} = $db_names[$i-1];
	    $ortho{$fam_id}{$sp_names[$i-1]}{$id} = 1;
	}
    }
}	

for $db( keys %gene_ids ){
    print STDERR "$db : ", scalar( keys %{$gene_ids{$db}} ), " entries\n";
}

## make a mapping of the sequence files to database identifiers
for $f(@seq_files){
    $sp = $f;
##    $sp =~ s/.+?\/(.+)$/$1/;
    $sp =~ s#.+?([^/]+)$#$1#;
    $sp =~ s/([^_]+[^\.]+).+$/$1/;
##    print STDERR "$f  --> $sp\n";
    $sp = lc($sp);
    if(!defined($db_map{$sp})){
	die "No database defined for species: $sp\n";
    }
    $seq_map{$sp} = $f;
}

## we now want to get the transcript information for each database..

## then we need to obtain for each gene in these lists
## 1. The canonical transcript coordinates
## 2. The genomic sequence covering the canonical transcript
## 3. The exons within that region, to enable splicing together of the canonical transcript
## 4. Write the transcript sequences to normal multiple sequence fasta files where each line
##    contains the sequence of one exon and each entry corresponds to an orthologue
##    Give some additional information in the id headers (although do not modify the normal id.

$gene_transcript_query = 
    "select a.stable_id as 'gene', c.stable_id as 'transcript', b.name as 'chr', a.seq_region_strand as 'strand', c.seq_region_start as 'start', c.seq_region_end as 'end', ".
    "e.exon_id as 'ex.id', e.seq_region_start as 'e.start', e.seq_region_end as 'e.end', d.rank ".
    "from gene a ".
    "inner join seq_region b on a.seq_region_id=b.seq_region_id ".
    "inner join transcript c on a.canonical_transcript_id=c.transcript_id ".
    "inner join exon_transcript d on c.transcript_id=d.transcript_id ".
    "inner join exon e on e.exon_id=d.exon_id ".
    "order by a.stable_id, d.rank;";

print STDERR "Getting exon data\n";
for $i(0..$#db_names){
    my $sp_name = $sp_names[$i];
    %{$gene_data{$sp_name}} = get_filtered_table($gene_transcript_query, $db_names[$i], 0, $gene_ids{$db_names[$i]});
    print STDERR "$db_names[$i] : ", scalar( @{$gene_data{$sp_name}{data}} ), "  rows\n";
    ## create a map of gene identifiers to row numbers in the given table.
    for $j(0..$#{$gene_data{$sp_name}{data}}){
	push @{$gene_exon{$sp_name}{ $gene_data{$sp_name}{data}[$j][0] }}, $j;
    }
}


## and at this point we have everything that we need apart from the sequences from the different genomes.
## We have to essentially read each genome in at a go and then extract the relevant sequences..
## This is a bit painful, but can be done..

## Let us now go through the orthologues and obtain the sequences.
## Unfortunately we cannot do this in the simple way as we cannot read all of the
## sequence data into memory. This means that we have to first create a set of
## files for the orthologues and then we need to go through and append those files with
## sequences from each genome.

## make the files..
for $orth( keys %ortho ){
    $orth_f_e = "orthologue_transcripts/exons_".$orth;
    $orth_f_i = "orthologue_transcripts/introns_".$orth;
    open(OUT, ">", $orth_f_e) || die "unable to create $orth_f $!\n";
    $orth_files_e{$orth} = $orth_f_e;
    close(OUT);
    open(OUT, ">", $orth_f_i) || die "unable to create $orth_f $!\n";
    $orth_files_i{$orth} = $orth_f_i;
    close(OUT);
}

open(my $stats, ">", "orthologue_transcripts/exon_intron_stats.csv") || die "unable to open exon_intron_stats.csv $!\n";

##@sp = keys %seq_map;
for $sp(keys %seq_map){
    print STDERR "getting sequences from genomes: $sp\n";
    ## let us make a map of the columns of the gene data set
    my %col_i = ();
    for my $j( 0..$#{$gene_data{$sp}{fields}} ){
	$col_i{ $gene_data{$sp}{fields}[$j] } = $j;
    }
    
    my %seq = read_fasta( $seq_map{$sp} );
    for $orth( keys %ortho ){
##	print STDERR "looking for $orth in $sp\n";
	next if(!defined( $ortho{$orth}{ $sp } ));
	my @gids = keys %{$ortho{$orth}{ $sp }};
##	print STDERR "Got : @gids\n";
	open(my $out, ">>", $orth_files_e{ $orth } ) || die "unable to open $orth_files_e{$orth} $!\n";
	open(my $out2, ">>", $orth_files_i{ $orth } ) || die "unable to open $orth_files_i{$orth} $!\n";
	for my $gid(@gids){
	    my @j = @{$gene_exon{$sp}{ $gid }};
	    my $header = ">".$gene_data{$sp}{data}[$j[0]][$col_i{'gene'}]."\t".$gene_data{$sp}{data}[$j[0]][$col_i{'transcript'}].
		"\t".$gene_data{$sp}{data}[$j[0]][ $col_i{'chr'} ]." ".$gene_data{$sp}{data}[$j[0]][$col_i{'start'}].
		"..".$gene_data{$sp}{data}[$j[0]][$col_i{'end'}].":".$gene_data{$sp}{data}[$j[0]][$col_i{'strand'}].
		"\t".$db_map{$sp}."\n";
	    print $out $header;
	    print $out2 $header;
	    print $stats $sp, "\t", $db_map{$sp}, "\t", $orth, "\t", $gene_data{$sp}{data}[$j[0]][$col_i{'gene'}], "\t",
		$gene_data{$sp}{data}[$j[0]][$col_i{'transcript'}], "\t", $gene_data{$sp}{data}[$j[0]][ $col_i{'chr'} ],
		"\t", $gene_data{$sp}{data}[$j[0]][$col_i{'strand'}], "\t", $gene_data{$sp}{data}[$j[0]][$col_i{'start'}], "..",
		$gene_data{$sp}{data}[$j[0]][$col_i{'end'}], "\t";
	    my $last_stop = -1;
	    my @e_sizes = ();
	    my @i_sizes = ();
	    for my $k( 0..$#j ){
		$j = $j[$k];
		## we should be able to directly get the exon information here..
		my $chr = $gene_data{$sp}{data}[$j][ $col_i{'chr'} ];
		my $start = $gene_data{$sp}{data}[$j][ $col_i{'e.start'} ];
		my $stop =  $gene_data{$sp}{data}[$j][ $col_i{'e.end'} ];
		my $strand = $gene_data{$sp}{data}[$j][ $col_i{'strand'} ];
		my $e_seq = substr($seq{$chr}, $start-1, 1 + $stop - $start);
		$e_seq = rev_comp($e_seq) if $strand == -1;
		push @e_sizes, length($e_seq);
		print $out $e_seq, "\n";
		if($last_stop >= 0){
		    my $i_seq = ($strand == 1) ? 
			substr($seq{$chr}, $last_stop, ($start - $last_stop)-1) :
			substr($seq{$chr}, $stop, ($last_stop - $stop)-1);
		    $i_seq = rev_comp($i_seq) if $strand == -1;
		    push @i_sizes, length($i_seq);
		    print $out2 $i_seq, "\n";
		}
		$last_stop = ($strand == 1) ? $stop : $start;
	    }
	    print $stats join(",", @e_sizes), "\t", join(",", @i_sizes), "\n";
	}
	close($out);
	close($out2);
    }
}

## then go through the sequence files and read in the genome sequences.. 

    
sub get_table {
    my($query, $db_name) = @_;
    my $dbs = "DBI:mysql:database=%s:host=ensembldb.ensembl.org:port=5306";
    my $db = DBI->connect(sprintf($dbs, $db_name), "anonymous", "") || die $DBI::errstr;
    my $sth = $db->prepare($query);
    $sth->execute();
    my $names = $sth->{NAME};
    my $numFields = $sth->{'NUM_OF_FIELDS'};
    my %data = ();
    $data{fields} = [@$names];
    my $line_no=0;
    while(my $ref = $sth->fetchrow_arrayref){
	for(my $i=0; $i < $numFields; $i++){
	    $data{data}[$line_no][$i] = defined($$ref[$i]) ? $$ref[$i] : "NULL";
	}
	$line_no++;
    }
    $db->disconnect();
    return(%data);
}

## returns only lines which have entries in the hash pointed to by id_set
sub get_filtered_table {
    my($query, $db_name, $col_i, $id_set) = @_;
##    my $dbs = "DBI:mysql:database=%s:host=ensembldb.ensembl.org:port=5306";
    my $dbs = "DBI:mysql:$db_name:mysql_read_default_file=$ENV{HOME}/.my.cnf";
    my $db = DBI->connect($dbs, "lmj") || die $DBI::errstr;
##    my $db = DBI->connect(sprintf($dbs, $db_name)) || die $DBI::errstr;
##    my $db = DBI->connect(sprintf($dbs, $db_name), "anonymous", "") || die $DBI::errstr;
    my $sth = $db->prepare($query);
    $sth->execute();
    my $names = $sth->{NAME};
    my $numFields = $sth->{'NUM_OF_FIELDS'};
    my %data = ();
    $data{fields} = [@$names];
    my $line_no=0;
    while(my $ref = $sth->fetchrow_arrayref){
	next if(!defined($id_set->{ $$ref[ $col_i ] }));
	for(my $i=0; $i < $numFields; $i++){
	    $data{data}[$line_no][$i] = defined($$ref[$i]) ? $$ref[$i] : "NULL";
	}
	$line_no++;
    }
    $db->disconnect();
    return(%data);
}

## read from a gzipped file
## id lines must end with REF
sub read_fasta {
    my $f = shift @_;
    my $in;
    if($f !~ /.gz$/){
	open($in, "<", $f) || die "unable to open $f $!\n";
    }else{
	open($in, "-|", "zcat $f") || die "unable to zcat $f\n";
    }
    my %seq = ();
    my $id = "";
    my $seq_is_ref = 0;
    while(<$in>){
	chomp;
	if($_ =~ /^>(\S+)/){
	    $id = $1;
	    $seq_is_ref = ($_ =~ /\sREF\s*$/);
	    if($seq_is_ref){
		$seq{$id} = "";
	    }
	    next;
	}
	if($seq_is_ref){
	    $seq{$id} .= $_;
	}
    }
    return(%seq);
}

sub rev_comp {
    my $seq = shift @_;
    $seq =~ tr/ACGTUMRWSYKVHDBN/TGCAAKYWSRMBDHVN/;
    $seq =~ tr/acgtumrwsykvhdbn/tgcaakywsrmbdhvn/;
    $seq = reverse($seq);
    return($seq);
}


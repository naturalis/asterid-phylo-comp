#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my ( $table, $outfile, $marker );
my $verbosity = WARN;
GetOptions(
	'table=s'   => \$table,
	'outfile=s' => \$outfile,
	'marker=s'  => \$marker,
	'verbose+'  => \$verbosity,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level'  => $verbosity,
	'-class'  => 'main',
);
my $io = Bio::SeqIO->new(
	'-format' => 'fasta',
	'-file'   => ">${outfile}",
);
my $gb = Bio::DB::GenBank->new;
$gb->retrieval_type('io_string');

# read the table
my %acc_by_spec;
{
	$log->info("going to read marker table $table");
	my ( $idx_spec, $idx_marker );
	open my $fh, '<', $table or die $!;
	LINE: while(<$fh>) {
		chomp;
		my @record = split /\t/, $_;
		if ( not defined $idx_spec and not defined $idx_marker ) {
			for my $i ( 0 .. $#record ) {
				$idx_spec   = $i if $record[$i] eq 'species';
				$idx_marker = $i if $record[$i] eq $marker;
			}
			next LINE;
		}
		my $species   = $record[$idx_spec];
		my $accession = $record[$idx_marker];
		if ( $accession ne 'NA' ) {
			$acc_by_spec{$species} = $accession;
		}
	}
	$log->info("read ".scalar(keys(%acc_by_spec))." accessions for $marker from $table");
}

# start fetching the sequences
$log->info("going to fetch sequences");
for my $species ( sort { $a cmp $b } keys %acc_by_spec ) {
	my $accession = $acc_by_spec{$species};
	my $defline = "${species}_${accession}";
	$defline =~ s/ /_/g;
	
	# do the fetch, write the seq
	my $seq = $gb->get_Seq_by_acc($accession);
	$seq->display_id($defline);
	$io->write_seq($seq);
	$log->info("wrote seq $defline");
}
$log->info("done");
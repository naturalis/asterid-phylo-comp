#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::AlignIO;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;

# process command line arguments
my ( $table, $outdir );
my $verbosity = WARN;
GetOptions(
	'table=s'  => \$table,
	'outdir=s' => \$outdir,
	'verbose+' => \$verbosity,
);

# instantiate helper objects
my $sgs = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# start reading the table
my ( @header );
$log->info("going to read marker table $table");
open my $fh, '<', $table or die $!;
LINE: while(<$fh>) {
	chomp;
	my @record = split /\t/;
	@header = @record and next LINE unless @header;
	
	# check if we have data
	if ( $record[1] eq 'NA' ) {
		$log->warn("no data for taxon $record[0]");
	}
	else {
		for my $i ( 3 .. $#record ) {
			my $m   = $header[$i];
			my $acc = $record[$i];
			if ( $acc ne 'NA' ) {
				$log->info("fetching $m sequence $acc for $record[1] ($record[0])");
				my $seq = $sgs->single_seq({ 'acc' => $acc });
				my $gi = $seq->gi;
				my $ti = $seq->ti;
			}
		}	
	}
}
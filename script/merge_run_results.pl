#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my ( $indir, $merged, $rerun, $limit, @markers );
GetOptions(
	'indir=s'   => \$indir,
	'super=s'   => \$merged,
	'rerun=s'   => \$rerun,
	'limit=i'   => \$limit,
	'markers=s' => \@markers,
	'verbose+'  => \$verbosity,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# read re-run table
my %rerun;
{
	$log->info("going to read rerun table $rerun");
	open my $fh, '<', $rerun or die $!;
	my @header;
	LINE: while(<$fh>) {
		chomp;
		my @record = split /\t/;
		@header = @record and next LINE unless @header;
		my $species = $record[1];
		my %record  = map { $header[$_] => $record[$_] } 1 .. $#header;
		$rerun{$species} = \%record;
	}
	$log->info("read ".scalar(keys(%rerun))." records");
}

# read merged table
my %merged;
{
	$log->info("going to read super table $merged");
	open my $fh, '<', $merged or die $!;
	my @header;
	LINE: while(<$fh>) {
		chomp;
		my @record = split /\t/;
		@header = @record and next LINE unless @header;
		my $taxon = $record[0];
		$merged{$taxon} = [] if not $merged{$taxon};
		
		# higher taxon had no data in initial run
		if ( $record[1] eq 'NA' ) {
			$log->info("$taxon had no data in initial run");
			
			# TNRS succeeded
			my $have_data;
			if ( -e "$indir/$taxon/species.tsv" ) {
				$log->info("going to read taxa table $indir/$taxon/species.tsv");
				open my $sfh, '<', "$indir/$taxon/species.tsv" or die $!;
				while(<$sfh>) {
					next if /^name/;
					my ($species) = split /\t/, $_;
					$species =~ s/_/ /g;
					$log->info("going to lookup species $species");	
					push @{ $merged{$taxon} }, $rerun{$species} if $rerun{$species};
					if ( $rerun{$species} ) {
						$log->info("adding record for $species");
						$have_data++;
					}
				}
			}
			
			# still no data after rerun
			if ( not $have_data ) {
				$log->warn("no data for $taxon after rerun");
				my %pseudo = map { $_ => "NA" } @markers;
				$pseudo{'species'} = "NA";
				$pseudo{'total'} = 0;
				push @{ $merged{$taxon} }, \%pseudo;			
			}		
		}
		else {
			my %record = map { $header[$_] => $record[$_] } 1 .. $#header;
			push @{ $merged{$taxon} }, \%record;
		}
		
		# recalculate counts with current marker set
		for my $r ( @{ $merged{$taxon} } ) {
			my $count = 0;
			$r->{$_} ne 'NA' and $count++ for @markers;
			$r->{'total'} = $count;
		}
	}
}

# print header
print join("\t",
	'taxa_in_original_tree',
	'species',
	'total',
	@markers
), "\n";

# print records
for my $t ( sort { $a cmp $b } keys %merged ) {
	my @sorted = sort { $b->{'total'} <=> $a->{'total'} } @{ $merged{$t} };
	my $i;
	RECORD: for my $r ( @sorted ) {
		print join("\t",
			$t,
			$r->{'species'},
			$r->{'total'},
			map { $r->{$_} } @markers
		), "\n";
		$i++;
		last RECORD if $limit and $i == $limit;
	}
}
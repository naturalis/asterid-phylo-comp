#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# process command line arguments
my ( $indir, @markers );
GetOptions(
	'indir=s'   => \$indir,
	'markers=s' => \@markers,
);

# start reading from directory
my %taxon;
opendir my $dh, $indir or die $!;
ENTRY: while( my $entry = readdir $dh ) {
	next ENTRY if $entry =~ /^\.\.?$/;
	next ENTRY unless -d "${indir}/${entry}";
	my @entries;
	
	# have a marker table
	if ( -e "${indir}/${entry}/markers-backbone.tsv" ) {
	
		# start reading marker table
		my @header;
		open my $fh, '<', "${indir}/${entry}/markers-backbone.tsv" or die $!;
		LINE: while(<$fh>) {
			chomp;
			my @record = split /\t/;
			
			# read the header if not yet done so
			@header = @record and next LINE unless @header;
			
			# read the record into a hash
			my %record = map { $header[$_] => $record[$_] } 0 .. $#header;
			
			# check if we have any of the markers
			my $flag;
			for my $m ( @markers ) {
				if ( $record{$m} ) {
					$flag++;
				}
				else {
					$record{$m} = "NA";
				}
			}
			$record{'total'} = $flag;
			push @entries, \%record if $flag;
		}
	}
	
	# sort by descending number of markers
	$taxon{$entry} = [ sort { $b->{'total'} <=> $a->{'total'} } @entries ];
	
	# no data found!
	if ( not @{ $taxon{$entry} } ) {
		my %pseudo = map { $_ => "NA" } @markers;
		$pseudo{'taxon'} = "NA";
		$pseudo{'total'} = 0;
		push @{ $taxon{$entry} }, \%pseudo;
	}
}

# print header
print join("\t",
	'taxa_in_original_tree', # higher taxon
	'species',               # species with at least one marker
	'total',                 # total markers present for species
	@markers                 # the list of markers
), "\n";

# print records
for my $t ( sort { $a cmp $b } keys %taxon ) {
	for my $record ( @{ $taxon{$t} } ) {
		print join("\t",
			$t,
			$record->{'taxon'},
			$record->{'total'},
			map { $record->{$_} } @markers,
		), "\n";
	}
}
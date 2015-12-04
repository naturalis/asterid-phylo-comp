#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my $taxon_id  = 'taxon_id'; # column in spreadsheet
my $original  = 'taxa_in_original_tree'; $ column in spreadsheet
my $species   = 'species'; # column in spreadsheet
my ( $intree, $data );
GetOptions(
	'intree=s'   => \$intree,
	'data=s'     => \$data,
	'verbose+'   => \$verbosity,
	'taxon_id=s' => \$taxon_id,
	'original=s' => \$original,
	'species=s'  => \$species,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level'  => $verbosity,
	'-class'  => 'main',
);
my $tree = parse_tree(
	'-format' => 'figtree',
	'-file'   => $intree,
);

# read species table
my %mapping;
{
	$log->info("going to read spreadsheet data from $data");
	my @header;
	my ( $s_idx, $t_idx, $c_idx ); # species, taxon, clade indices
	open my $fh, '<', $data or die $!;
	LINE: while(<$fh>) {
		chomp;
		my @line = split /\t/, $_;
		if ( not @header ) {
			@header = @line;
			for my $i ( 0 .. $#header ) {
				$s_idx = $i if $header[$i] eq $species;
				$c_idx = $i if $header[$i] eq $original;
				$t_idx = $i if $header[$i] eq $taxon_id;
			}
			next LINE;
		}
		$mapping{ $line[$t_idx] } = [ $line[$c_idx], $line[$s_idx] ];	
	}
}

# remap tree
$tree->visit_dept_first(
	'-post' => sub {
		my $n = shift;
		my @c = @{ $n->get_children };
		
		# node is a tip
		if ( not @c ) {
			my $name = $n->get_name;
			if ( $mapping{$name} ) {
				$n->set_name( $mapping{$name}->[1] );
				$n->set_generic( 'clade' => $mapping{$name}->[0] );				
			}
			else {
				$log->error("no mapping for $name");
			}
		}
		
		# node is internal, carry over the clade grouping
		else {
		
		}	
	}
);
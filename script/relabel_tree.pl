#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util 'sum';
use Bio::Phylo::Treedrawer;
use Bio::Phylo::IO qw'parse unparse';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::CONSTANT ':objecttypes';

# process command line arguments
my $verbosity = WARN;
my $taxon_id  = 'taxon_id'; # column in spreadsheet
my $original  = 'taxa_in_original_tree'; # column in spreadsheet
my $species   = 'species'; # column in spreadsheet
my $character = 'perforation plate type'; # column in spreadsheet
my $posterior = 0.5; # threshold posterior
my ( $intree, $data, $ancstates );
GetOptions(
	'intree=s'    => \$intree,
	'data=s'      => \$data,
	'verbose+'    => \$verbosity,
	'taxon_id=s'  => \$taxon_id,
	'original=s'  => \$original,
	'species=s'   => \$species,
	'posterior=f' => \$posterior,
	'ancstates=s' => \$ancstates,
	'character=s' => \$character,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level'       => $verbosity,
	'-class'       => 'main',
);
my $project = parse(
	'-format'      => 'figtree',
	'-file'        => $intree,
	'-as_project'  => 1,
);
my ($tree) = @{ $project->get_items(_TREE_) };
my $drawer = Bio::Phylo::Treedrawer->new(
	'-format'      => 'svg',
	'-mode'        => 'clado',
	'-shape'       => 'rect',
	'-width'       => 1600,
	'-height'      => 4000,
	'-tree'        => $tree,
	'-node_radius' => 10,
	'-tip_radius'  => 10,
	'-text_width'  => 600,
	'-text_horiz_offset' => 15,
);

# read species table
my ( %mapping, %char );
{
	$log->info("going to read taxon data from $data");
	my @header;
	my ( $s_idx, $t_idx, $c_idx, $p_idx ); # species, taxon, clade, phenotype indices
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
				$p_idx = $i if $header[$i] eq $character;
			}
			next LINE;
		}
		$mapping{ $line[$t_idx] } = "'" . $line[$s_idx] . ' / ' . $line[$c_idx] . "'";	
		$char{ $line[$t_idx] } = $line[$p_idx];
	}
}

# apply ancestral state pies
{
	
	# read in the state reconstructions
	$log->info("going to read ancestral states from $ancstates");
	my ( %header, %record );
	open my $fh, '<', $ancstates or die $!;
	LINE: while(<$fh>) {
		chomp;
		my @line = split /\t/, $_;
		if ( not %header ) {
			for my $i ( 0 .. $#line ) {
				if ( $line[$i] =~ /Node(\d+) P\((\d)\)/ ) {
					my ( $index, $state ) = ( $1, $2 );
					$header{$index} = {} if not $header{$index};
					$header{$index}->{$state} = $i;
					$record{$index} = {} if not $record{$index};
					$record{$index}->{$state} = [];
				}
			}
			next LINE;
		}
		for my $index ( keys %record ) {
			for my $state ( keys %{ $record{$index} } ) {
				my $i = $header{$index}->{$state};
				if ( $line[$i] ne '--' ) {
					push @{ $record{$index}->{$state} }, $line[$i];
				}
			}
		}	
	}
	
	# apply state reconstructions to the generic 'pie' slot
	my $index = 1;
	$tree->visit_depth_first(
		'-post' => sub {
			my $n = shift;
			$n->set_font_face('Verdana');
			if ( $n->is_internal ) {
				my $values = $record{$index};
				for my $key ( keys %{ $values } ) {
					my @v = @{ $values->{$key} };
					$values->{$key} = sum(@v) / scalar(@v);					
				}
				$n->set_generic( 'pie' => $values );
				$n->set_node_outline_color( 'black' );
				$index++;
			}
			else {
				$n->set_font_style('Italic');
			}
		}
	);
}

# synonyms that slipped into @fredericlens's spreadsheet
my %synonym = (
	'1317872' => '59675',
	'1630344' => '1630342',
	'237952'  => '85293',
);

# relabel tips
my %colors;
for my $tip ( @{ $tree->get_terminals } ) {
	my $name = $tip->get_name;
	my $key = $mapping{$name} ? $name : $synonym{$name};
	if ( $mapping{$key} ) {
		$tip->set_name($mapping{$key});
		my $state = $char{$key};
		if ( not %colors ) {
			$colors{$state} = 'white';
		}
		elsif ( not $colors{$state} ) {
			$colors{$state} = 'black';
		}
		$tip->set_node_color($colors{$state});
		$tip->set_node_outline_color( 'black' );
	}
	else {
		$log->error("no mapping for $name");
	}
}

# collapse low support nodes
$tree->visit_depth_first(
	'-post' => sub {
		my $n = shift;
		$n->set_font_size('12px');
		if ( $n->is_internal ) {
			my $p = $n->get_meta_object('fig:posterior');
			if ( $p < $posterior ) {
				$n->collapse;
			}
			else {
				$n->set_name( sprintf "%.2f", $p );
			}
		}
	}
);
$tree->ladderize;

print $drawer->draw;
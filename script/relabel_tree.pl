#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util 'sum';
use List::MoreUtils 'uniq';
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
my @excise    = ( 
	'Poraqueiba guianensis',
	'Alstonia scholaris',
	'Maesa indica',
	'Berzelia lanuginosa',
	'Brunia albiflora',
	'Abelia spathulata',
);
my ( $intree, $data, $ancstates );
GetOptions(
	'intree=s'     => \$intree,
	'data=s'       => \$data,
	'verbose+'     => \$verbosity,
	'taxon_id=s'   => \$taxon_id,
	'original=s'   => \$original,
	'species=s'    => \$species,
	'posterior=f'  => \$posterior,
	'ancstates=s'  => \$ancstates,
	'character=s'  => \$character,
	'excise=s'     => \@excise,
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
	'-height'      => 2000,
	'-tree'        => $tree,
	'-node_radius' => 10,
	'-tip_radius'  => 10,
	'-text_width'  => 600,
	'-text_horiz_offset'     => 15,
	'-collapsed_clade_width' => 2,
);

# main operations
my %meta   = read_species_table($data);
my %record = read_ancestral_states($ancstates);
my %colors = apply_states( \%record, \%meta );
collapse_tree( \@excise, \%colors );

# read species table
sub read_species_table {
	my $data = shift;

	$log->info("going to read taxon data from $data");
	my @header;
	my %meta;
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
		$meta{ $line[$t_idx] } = {
			'species' => $line[$s_idx],
			'clade'   => $line[$c_idx],
			'state'   => $line[$p_idx],
		};
	}
	return %meta;
}

# read in the state reconstructions
sub read_ancestral_states {
	my $ancstates = shift;
	
	$log->info("going to read ancestral states from $ancstates");
	my ( %header, %record );
	open my $fh, '<', $ancstates or die $!;
	LINE: while(<$fh>) {
		chomp;
		my @line = split /\t/, $_;
		
		# parse the header, instantiate record fields
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
		
		# store the record
		for my $index ( keys %record ) {
			for my $state ( keys %{ $record{$index} } ) {
				my $i = $header{$index}->{$state};
				if ( $line[$i] ne '--' ) {
					push @{ $record{$index}->{$state} }, $line[$i];
				}
			}
		}	
	}
	return %record;
}

# apply ancestral state pies
sub apply_states {
	my ( $record, $meta ) = @_;
	$log->info("going to apply states");

	# synonyms that slipped into @fredericlens's spreadsheet
	my %synonym = (
		'1317872' => '59675',
		'1630344' => '1630342',
		'237952'  => '85293',
	);

	# apply state reconstructions to the generic 'pie' slot
	my $index = 1;
	my %colors;
	$tree->visit_depth_first(
		'-post' => sub {
			my $n = shift;
			
			# both tips and nodes
			$n->set_font_face('Verdana');
			$n->set_font_size('12px');
			
			# nodes have multiple states in 'pie' slot
			if ( $n->is_internal ) {
				my $values = $record->{$index};
				
				# compute averages for each state
				for my $key ( keys %{ $values } ) {
					my @v = @{ $values->{$key} };
					$values->{$key} = sum(@v) / scalar(@v);					
				}
				$n->set_generic( 'pie' => $values );
				$n->set_node_outline_color( 'black' );
				$index++;
			}
			else {
			
				# tips are species names, hence italic
				$n->set_font_style('Italic');
				my $name = $n->get_name;
				my $key = $meta->{$name} ? $name : $synonym{$name};		
				if ( $meta->{$key} ) {
					my $state = $meta->{$key}->{'state'};
					
					# progressively assign colors
					if ( not %colors ) {
						$colors{$state} = 'white';
					}
					elsif ( not $colors{$state} ) {
						$colors{$state} = 'black';
					}	
					
					# set name / clade / state				
					$n->set_name( $meta->{$key}->{'species'} );
					$n->set_generic( 'clade' => [ $meta->{$key}->{'clade'} ] );
					$n->set_generic( 'state' => [ $state ] );
					$n->set_node_color($colors{$state});
					$n->set_node_outline_color( 'black' );					
				}
				else {
					$log->error("no mapping for $name");
				}		
			}
		}
	);
	return %colors;		
}

# collapse low-support nodes, unwanted tips and monophyletic,
# uniform-state clades
sub collapse_tree {
	my ( $excise, $colors ) = @_;
	$log->info("going to collapse tree");

	# collapse low support nodes
	$tree->visit_depth_first(
		'-post' => sub {
			my $n = shift;		
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

	# prune unwanted tips
	my @prune;
	for my $name ( @$excise ) {
		if ( my $tip = $tree->get_by_name($name) ) {
			push @prune, $tip;
		}
		else {
			$log->error("Can't find $name in tree");
		}
	}
	$tree->prune_tips(\@prune);

	# collapse monophyletic, uniform-state clades
	$tree->visit_depth_first(
		'-post' => sub {
			my $n = shift;
			my @c = @{ $n->get_children };
			if ( @c ) {
				my ( @clade, @state, $size );
				for my $c ( @c ) {
					push @clade, @{ $c->get_generic('clade') };
					push @state, @{ $c->get_generic('state') };
					$size += $c->get_generic('size') || 1;
				}
				$n->set_generic( 'clade' => \@clade );
				$n->set_generic( 'state' => \@state );
				$n->set_generic( 'size'  => $size );
				if ( scalar(uniq(@clade)) == 1 and scalar(uniq(@state)) == 1 ) {
					my ($clade) = @clade;
					my ($state) = @state;
					my $prob    = sprintf '%.2f', $n->get_meta_object('fig:posterior');
					my $name = $clade . " (n=$size, p=$prob)";
					$n->set_collapsed(1);
					$n->set_name($name);
					$n->set_node_color($colors->{$state});
				}
			}
		}
	);

	$tree->ladderize;
}

$log->info("going to draw tree");
print $drawer->draw;
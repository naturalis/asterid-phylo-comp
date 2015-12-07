#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';

# process command line arguments
my $burnin  = 0.1;
my $outdir  = '.';
my $species = 'taxon_id'; # column in final_selection_table_asterids.tsv
my $char    = 'perforation plate type'; # column in final_selection_table_asterids.tsv
my ( $tree, $data, $consensus );
GetOptions(
	'tree=s'      => \$tree,
	'data=s'      => \$data,
	'burnin=f'    => \$burnin,
	'outdir=s'    => \$outdir,
	'species=s'   => \$species,
	'consensus=s' => \$consensus,
);

my %synonym = (
	'59675'   => '1317872',
	'1630342' => '1630344',
	'85293'   => '237952',
);

# read all trees into array
my @lines;
{
	open my $fh, '<', $tree or die $!;
	while(<$fh>) {
		chomp;
		push @lines, $_ if /\S/;
	}
}

# subset the array to discard the burnin
my $start = int( $burnin * scalar @lines );
my @trees = @lines[$start .. $#lines];

# create the translation table
my %translate;
my $index = 1;
my ($first) = @trees;
{
	my $t = parse_tree( '-format' => 'newick', '-string' => $first );
	for my $tip ( @{ $t->get_terminals } ) {
		my $name = $tip->get_name;
		$translate{$name} = $index++;
	}
}

# translate the tip labels according to the translation table
for my $t ( @trees ) {
	while( my ( $key, $value ) = each %translate ) {
		$t =~ s/([,\(])$key:/$1$value:/;
	}
}

# write the tree file
{
	my $outfile = $tree;
	$outfile =~ s/.+\///;
	$outfile =~ s/\.[^\.]+//;
	$outfile = $outdir . '/' . $outfile . '.trees';
	open my $outfh, '>', $outfile or die $!;
	
	# nexus header
	print $outfh '#NEXUS', "\n", 'BEGIN TREES;', "\n", "\tTRANSLATE\n";
	
	# translation table
	my $ntax = scalar keys %translate;
	for my $label ( sort { $translate{$a} <=> $translate{$b} } keys %translate ) {
		my $index  = $translate{$label};
		my $suffix = $index == $ntax ? "\n\t\t;" : ',';
		print $outfh "\t\t", $index, "\t", $label, $suffix, "\n";	
	}
	
	# trees
	my $index = 1;
	for my $t ( @trees ) {
		print $outfh "\tTREE TREE_$index = $t\n";
		$index++;
	}
	
	# nexus footer
	print $outfh 'END;';
}

# read and write the data file
{
	open my $fh, '<', $data or die $!;
	my $outfile = $data;
	$outfile =~ s/.+\///;
	$outfile =~ s/\.[^\.]+//;	
	$outfile = $outdir . '/' . $outfile . '.tsv';
	open my $outfh, '>', $outfile or die $!;
	
	# start reading the infile
	my @header;
	my ( $s_idx, $c_idx, %char );
	my $state = 0;
	LINE: while(<$fh>) {
		chomp;
		my @line = split /\t/, $_;
		
		# read the TSV header, find column indices for species and character
		if ( not @header ) {
			@header = @line;
			for my $i ( 0 .. $#header ) {
				$s_idx = $i if $header[$i] eq $species;
				$c_idx = $i if $header[$i] eq $char;
			}
			next LINE;
		}
		
		# only print records that occur in the tree
		my ( $taxon, $char ) = ( $line[$s_idx], $line[$c_idx] );
		$char = $char{$char} // ( $char{$char} = $state++ );
		if ( $translate{$line[$s_idx]} ) {
			print $outfh $taxon, "\t", $char, "\n";
		}
		else {
			if ( $synonym{$taxon} and $translate{$synonym{$taxon}} ) {
				print $outfh $synonym{$taxon}, "\t", $char, "\n";
			}
			else {
				die $taxon;
			}
		}
	}
}

# write the command file
{

	# open handle, write first lines
	open my $outfh, '>', $outdir . '/' . 'commands.txt' or die $!;
	print $outfh "1\n"; # multistate
	print $outfh "2\n"; # MCMC
	print $outfh "PriorAll exp 10\n"; # priors for states
	
	# read the consensus tree, will estimate ancestral states for nodes in it
	my $t = parse_tree(
		'-format' => 'figtree',
		'-file'   => $consensus,
	);
	
	# write node statements
	my $index = 1;
	$t->visit_depth_first(
		'-post' => sub {
			my $n = shift;
			my @c = @{ $n->get_children };
			if ( @c ) {
				my %tips;
				for my $c ( @c ) {
					my @t = @{ $c->get_generic('tips') };
					$tips{$_}++ for @t;
				}
				my @tips = keys %tips;
				print $outfh "AddNode Node" . $index++ . " @tips\n";
				$n->set_generic( 'tips' => [ @tips ] );
			}
			else {
				my $name = $n->get_name;
				my $ti = $translate{$name} || $translate{$synonym{$name}};
				$n->set_generic( 'tips' => [ $ti ] );
			}
		}
	);
	print $outfh "Run\n";
}
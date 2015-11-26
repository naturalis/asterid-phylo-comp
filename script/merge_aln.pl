#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::IO qw'parse_matrix unparse';

# process command line arguments
my ( @infile, $format );
GetOptions(
	'infile=s' => \@infile,
	'format=s' => \$format,
);

# instantiate helper objects
my $fac = Bio::Phylo::Factory->new;
my $project = $fac->create_project;
my $taxa    = $fac->create_taxa;
my $merged  = $fac->create_matrix( '-type' => 'dna', '-taxa' => $taxa );
$project->insert($taxa);
$project->insert($merged);

# read matrices, collect distinct taxa
my %taxa;
my @matrices;
for my $i ( @infile ) {
	my $matrix = parse_matrix(
		'-format' => 'fasta',
		'-type'   => 'dna',
		'-file'   => $i,
	);
	
	# clean up row names: only using binomials
	$matrix->visit(sub{
		my $row  = shift;
		my $name = $row->get_name;
		$name =~ s/ .+//;
		$name =~ s/_[^_]+$//;
		$row->set_name($name);
		$taxa{$name}++;
	});
	push @matrices, $matrix;
}

# populate the merged matrix
for my $name ( sort { $a cmp $b } keys %taxa ) {

	# assemble the concatenated sequence
	my $seq;
	for my $m ( @matrices ) {
		if ( my $row = $m->get_by_name($name) ) {
			my $char = $row->get_char;
			$seq .= $char;
		}
		else {
			my $nchar = $m->get_nchar;
			$seq .= '?' x $nchar;
		}
	}
	
	# create the taxon and datum objects
	my $taxon = $fac->create_taxon( '-name' => $name );
	$taxa->insert($taxon);	
	$merged->insert($fac->create_datum( 
		'-type'  => 'dna',
		'-name'  => $name,
		'-taxon' => $taxon,
		'-char'  => $seq,
	));
}

print unparse(
	'-phylo'  => $project,
	'-format' => $format,
);
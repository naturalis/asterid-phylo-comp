#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::GenBank;
use Bio::Phylo::Factory;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::IO qw'parse_matrix unparse';

# process command line arguments
my $format = 'phylip';
my $verbosity = WARN;
my @infile;
GetOptions(
	'infile=s' => \@infile,
	'format=s' => \$format,
	'verbose+' => \$verbosity,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new( '-level' => $verbosity, '-class' => 'main' );
my $dbg = Bio::DB::GenBank->new;
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
		my $row = shift;
		my $accession = $row->get_name;
		$accession =~ s/ .+//; # strip trailing defline
		$accession =~ s/^.+_([^_]+)$/$1/; # keep accession number
		
		# attempt to fetch taxon ID from sequence
		$log->info("attempting to fetch sequence $accession");
		my $seq = $dbg->get_Seq_by_acc($accession);
		my $tid = $seq->species->ncbi_taxid;
		$row->set_name($tid);
		$taxa{$tid}++;
	});
	push @matrices, $matrix;
}

# populate the merged matrix
NAME: for my $name ( sort { $a cmp $b } keys %taxa ) {

	# assemble the concatenated sequence, insert missing if need be
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
	
	# create the taxon and datum objects, insert in containers
	my $taxon = $fac->create_taxon( '-name' => $name );
	$taxa->insert($taxon);	
	$merged->insert($fac->create_datum( 
		'-type'  => 'dna',
		'-name'  => $name,
		'-taxon' => $taxon,
		'-char'  => $seq,
	));
}

# write output, can be nexus / phylip / nexml / fasta
print unparse(
	'-phylo'  => $project,
	'-format' => $format,
);
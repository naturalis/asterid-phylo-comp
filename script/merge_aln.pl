#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
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
my %ti_for_name;
my @matrices;
for my $i ( @infile ) {
	$log->info("going to read $i");
	my $matrix = parse_matrix(
		'-format' => 'fasta',
		'-type'   => 'dna',
		'-file'   => $i,
	);
	
	# parse out accession numbers, set as row name, 
	# store mapping to binomial
	my %name_for_accession;
	$matrix->visit(sub{
		my $row = shift;
		my $accession = $row->get_name;
		$accession =~ s/ .+//; # strip trailing defline
		my @parts = split /_/, $accession;
		$accession = pop @parts;
		my $name = join ' ', @parts;
		$name_for_accession{$accession} = $name;
		$log->info("$name -> $accession");
		$row->set_name($accession);
	});
	
	# stream all results	
	my $io = $dbg->get_Stream_by_acc([ keys %name_for_accession ]);
	while( my $seq = $io->next_seq ) {
		my $accession = $seq->accession_number;
		my $name = $name_for_accession{$accession};
		my $tid  = $ti_for_name{$name} || ($ti_for_name{$name} = $seq->species->ncbi_taxid);
		$log->info("$name -> $accession -> $tid");
		$taxa{$tid}++;
		my $row = $matrix->get_by_name($accession);
		$row->set_name($tid);
	}	
	push @matrices, $matrix;
}
$log->info("\n".Dumper(\%taxa)."\n");

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
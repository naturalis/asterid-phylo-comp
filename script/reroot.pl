#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my ( $infile, $outgroup );
my $verbosity = WARN;
GetOptions(
	'infile=s'   => \$infile,
	'outgroup=s' => \$outgroup,
	'verbose+'   => \$verbosity,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);
my @og = split /,/, $outgroup;

# read from infile
my $i = 1;
open my $fh, '<', $infile or die $!;
while(<$fh>) {
	chomp;
	my $newick = $_;
	my $tree = parse_tree(
		'-format' => 'newick',
		'-string' => $newick,
	);
	my @tips = map { $tree->get_by_name($_) } @og;
	my $tc = scalar(@tips);	
	my $mrca = $tree->get_mrca(\@tips);
	my $cc = scalar(@{$mrca->get_terminals});	
	if ( $cc != $tc ) {
		$log->warn("outgroup not monophyletic ($cc != $tc) in tree $i: ".$mrca->to_newick);
	}
	$mrca->set_root_below(100);
	print $tree->to_newick, "\n";
	$i++;
}
	
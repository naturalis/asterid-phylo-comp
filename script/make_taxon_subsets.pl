#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# process command line arguments
my $infile;
my $outdir;
GetOptions(
	'infile=s'   => \$infile,
	'outdir=s'   => \$outdir,
);

# read TSV
{
	my @header;
	open my $fh, '<', $infile or die $!;
	LINE: while(<$fh>) {
		chomp;
		my @record = split /\t/, $_;
		if ( not @header ) {
			@header = @record;
			next LINE;
		}
		if ( $record[1] ) {
			$record[2] =~ s/"//g;
			my @subtaxa = grep { /\S/ } split /,\s*/, $record[2];
			for my $taxon ( @subtaxa ) {
				if ( not -d "${outdir}/${taxon}" ) {
					mkdir "${outdir}/${taxon}";
				}
			}
		}
	}
}
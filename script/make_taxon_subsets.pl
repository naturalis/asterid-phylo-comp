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
			my $root = $record[0];
			$root =~ s/ /_/g;
			$root =~ s/"//g;
			if ( not -d "${outdir}/${root}" ) {
				mkdir "${outdir}/${root}";
			}
			if ( not -e "${outdir}/${root}/names.txt" ) {
				my @subtaxa = grep { /\S/ } split /,\s*/, $record[2];
				open my $fh, '>', "${outdir}/${root}/names.txt" or die $!;
				print $fh join ",", @subtaxa;
			}
		}
	}
}
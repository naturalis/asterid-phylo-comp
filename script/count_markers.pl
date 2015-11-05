#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# process command line arguments
my $indir;
GetOptions(
	'indir=s' => \$indir,
);

my %marker;
opendir my $dh, $indir or die $!;
while( my $entry = readdir $dh ) {
	next if $entry =~ /\.\.?/;
	
	# directory has a marker table
	if ( -d "${indir}/${entry}" && -e "${indir}/${entry}/markers-backbone.tsv" ) {
	
		# start reading the table
		my @header;
		open my $fh, '<', "${indir}/${entry}/markers-backbone.tsv" or die $!;
		LINE: while(<$fh>) {
			chomp;
			my @record = split /\t/;
			@header = @record and next LINE unless @header;
			for ( 1 .. $#header ) {
				$marker{$header[$_]}++ if $record[$_];
			}
		
		}
	}
}
print "marker\ttaxa\n";
print $_, "\t", $marker{$_}, "\n" for sort { $marker{$b} <=> $marker{$a} } keys %marker;
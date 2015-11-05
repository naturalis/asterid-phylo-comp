#!/bin/bash

SPREADSHEET=data/taxa/sampling_asterids_tree.tsv
SUPERSMART=data/supersmart

# make top level working directory
if [ ! -d "$SUPERSMART" ]; then
	mkdir $SUPERSMART
fi

# make sub directories
perl script/make_taxon_subsets.pl -i $SPREADSHEET -o $SUPERSMART

# make list of taxon names
TAXA=`ls $SUPERSMART/`

# iterate over taxa
for TAXON in $TAXA; do

	# move into working directory
	cd $SUPERSMART/$TAXON

	# perform taxonomic name resolution
	if [ ! -e "species.tsv" ]; then
		smrt taxize -b -r $TAXON
	fi
	
	# align all phylota clusters for the species
	if [ ! -e "aligned.txt" ]; then
		smrt align
	fi
	
	# run all-vs-all BLAST merger on phylota clusters
	if [ ! -e "merged.txt "]; then
		smrt orthologize
	fi
	
	# perform a supermatrix merger, produce marker table
	if [ ! -e "markers-backbone.tsv" ]; then
		smrt bbmerge --exemplars -1
	fi
	
	# move back
	cd -
done

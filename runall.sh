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

	# run smrt taxize
	if [ ! -e "species.tsv" ]; then
		smrt taxize -b -r $TAXON
	fi
	
	# run smrt align
	if [ ! -e "aligned.txt" ]; then
		smrt align
	fi
	
	# move back
	cd -

done
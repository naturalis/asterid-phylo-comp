#!/bin/bash

SPREADSHEET=data/taxa/sampling_asterids_tree.tsv
SUPERSMART=data/supersmart-merged
MARKERS=data/markers
export SUPERSMART_BACKBONE_MAX_DISTANCE=0.25
export SUPERSMART_BACKBONE_MIN_COVERAGE=1

# make top level working directory
if [ ! -d "$SUPERSMART" ]; then
	mkdir $SUPERSMART
fi

# make sub directories
perl script/make_taxon_subsets.pl -i $SPREADSHEET -o $SUPERSMART

# make list of taxon names
TAXA=`ls $SUPERSMART/`

# iterate over taxa to create marker tables
for TAXON in $TAXA; do

	# move into working directory
	cd $SUPERSMART/$TAXON

	# perform taxonomic name resolution
	if [ ! -e "species.tsv" ]; then
		smrt taxize -b -r `cat names.txt`
	fi
	
	# align all phylota clusters for the species
	if [ ! -e "aligned.txt" ]; then
		smrt align				
	fi	
	
	# run all-vs-all BLAST merger on phylota clusters
	if [ -e "aligned.txt" ] && [ ! -e "merged.txt" ]; then
		smrt orthologize
	else
		if [ -e "merged.txt" ]; then
			echo "already clustered $TAXON"
		else
			echo "no alignments were made for $TAXON, won't cluster"
		fi
	fi	
	
	# perform a supermatrix merger, produce marker table
	if [ -e "merged.txt" ] && [ ! -e "markers-backbone.tsv" ]; then
		smrt bbmerge --exemplars -1
	else
		if [ -e "markers-backbone.tsv" ]; then
			echo "already merged $TAXON"
		else
			echo "no clustering was done for $TAXON, won't merge"
		fi
	fi	
	
	# move back
	cd -
done

# create marker summary
if [ ! -e "$MARKERS/markers.tsv" ]; then
	perl script/count_markers.pl -i $SUPERSMART > "$MARKERS/markers.tsv"
fi
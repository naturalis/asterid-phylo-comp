#!/bin/bash

SPREADSHEET=data/taxa/sampling_asterids_tree.tsv
SUPERSMART=data/supersmart-merged
MARKERS=data/markers
ENRICHED=data/enriched
ENRICHED_MARKERS=$ENRICHED/final_selection_table_asterids.txt
RERUN=data/supersmart-rerun
ALIGNMENTS=data/alignments
SCRIPT=`pwd`/script
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


##########################################################################################
# Here we do an initial run of the SUPERSMART pipeline where we try to collect data for
# each of the higher taxa in the input spreadsheet individually. This sometimes fails,
# especially in cases where an input taxon has <3 species in our database, i.e. the 
# pipeline considers such an input taxon near-monotypic, such that no phylogenetically
# informative data set can be constructed. In a second step, we will create a set of names
# only with species from these input taxa, so that we can then convince the pipeline that,
# in combination, phylogenetically informative data sets can still be constructed.

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

##########################################################################################
# Here we summarize the results of the initial run into a table with candidate markers,
# sorted by number of species for which we have data for this marker, and a table with
# species and a more or less arbitrary selection of frequently-sequenced markers.

# create marker summary
if [ ! -e "$MARKERS/markers.tsv" ]; then
	perl $SCRIPT/count_markers.pl -i $SUPERSMART > "$MARKERS/markers.tsv"
fi

# create taxon summary
if [ ! -e "$MARKERS/taxa.tsv" ]; then
	perl $SCRIPT/count_taxa.pl -i $SUPERSMART -m 'rbcL' -m 'matK' -m 'rps16' -m 'ndhF' \
	-m '5.8S-ribosomal-RNA' -m 'contains-trnT-trnL-intergenic-spacer,-tRNA-Leu-(trnL),-trnL-trnF-intergenic-spacer,-and-tRNA-Phe-(trnF)' \
	> "$MARKERS/taxa.tsv"
fi

##########################################################################################
# Now we can assess which taxa failed to generate data sets in the initial run. We first
# combine the species spreadsheets for these higher taxa, which we then use as input to
# re-run the pipeline, whose results we then summarize in the same way as above.

# create species table of missing taxa
if [ ! -e "$RERUN/species.tsv" ]; then
	MISSING=`awk '{if ($2 == "NA") print $1}' "$MARKERS/taxa.tsv"`
	for TAXON in $MISSING; do
		if [ ! -e "$RERUN/species.tsv" ]; then
			head -1 "$SUPERSMART/$TAXON/species.tsv" > "$RERUN/species.tsv"
		fi	
		grep -v '^name' "$SUPERSMART/$TAXON/species.tsv" >> "$RERUN/species.tsv"
	done
fi

# re-run SUPERSMART
cd $RERUN

# align all phylota clusters for the species
if [ ! -e "aligned.txt" ]; then
	smrt align				
fi	

# run all-vs-all BLAST merger on phylota clusters
if [ -e "aligned.txt" ] && [ ! -e "merged.txt" ]; then
	smrt orthologize
else
	if [ -e "merged.txt" ]; then
		echo "already clustered $RERUN"
	else
		echo "no alignments were made for $RERUN, won't cluster"
	fi
fi	

# perform a supermatrix merger, produce marker table
if [ -e "merged.txt" ] && [ ! -e "markers-backbone.tsv" ]; then
	smrt bbmerge --exemplars -1
else
	if [ -e "markers-backbone.tsv" ]; then
		echo "already merged $RERUN"
	else
		echo "no clustering was done for $RERUN, won't merge"
	fi
fi	

# create taxon summary
if [ ! -e "taxa.tsv" ]; then
	perl $SCRIPT/count_taxa.pl -i "markers-backbone.tsv" -m 'rbcL' -m 'matK' -m 'rps16' -m 'ndhF' \
	-m '5.8S-ribosomal-RNA' -m 'contains-trnT-trnL-intergenic-spacer,-tRNA-Leu-(trnL),-trnL-trnF-intergenic-spacer,-and-tRNA-Phe-(trnF)' \
	> "taxa.tsv"
fi

# move back
cd -

# merge the two run results
if [ ! -e "$MARKERS/merged.tsv" ]; then
	perl $SCRIPT/merge_run_results.pl -i $SUPERSMART -s $MARKERS/taxa.tsv -r $RERUN/taxa.tsv \
	-limit 3 -m 'rbcL' -m 'matK' -m 'rps16' -m 'ndhF' > $MARKERS/merged.tsv
fi

##########################################################################################
# Here we first fetch all the sequences for the selected markers and taxa in separate 
# files. We then align these and concatenate them as input for exabayes, i.e. in a PHYLIP
# file with short taxon names
MARKER_SELECTION="rbcL matK rps16 ndhF"
for M in $MARKER_SELECTION; do

	# do the download
	if [ ! -e "$ENRICHED/$M.fa" ]; then
		perl script/fetch_seqs.pl -t $ENRICHED_MARKERS -m $M -o $ENRICHED/$M.fa -v
	fi
	
	# do the alignment
	if [ ! -e "$ENRICHED/$M.aln.fa" ]; then
		mafft --auto $ENRICHED/$M.fa > $ENRICHED/$M.aln.fa
	fi
done

# do the concatenation
if [ ! -e "$ENRICHED/merged.phy" ]; then
	perl script/merge_aln.pl -i "$ENRICHED/rbcL.aln.fa" -i "$ENRICHED/matK.aln.fa" \
		"$ENRICHED/rps16.aln.fa" -i "$ENRICHED/ndhF.aln.fa" -f phylip > "$ENRICHED/merged.phy"
fi
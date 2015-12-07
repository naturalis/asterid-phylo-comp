# asterid-phylo-comp
Project repository to infer bayesian phylogenies of Asterids and to perform 
comparative analysis of state transitions in vascular morphology.

The following sections outline the analysis steps taken.

## Initial data mining
To obtain a first impression of which markers are, in combination, sufficiently 
sampled to cover the major groups of interest and, conversely, to identify which
taxa might be suitable representatives for their respective major groups we ran
the following steps of the [SUPERSMART](http://www.supersmart-project.org) pipeline:

1. `smrt taxize` as applied to each of the higher taxon names in the `names_to_query`
   column of the [sampling_asterids_tree](data/taxa/sampling_asterids_tree.tsv) spreadsheet,
   binned by their grouping under `taxa_in_original_tree`. For example, for the family
   Fouquieriaceae we expanded the genus Fouquieria.
2. For each of these `taxa_in_original_tree` we then ran `smrt align` and `smrt orthologize`
   to obtain their available PhyLoTA clusters, and then ran `smrt bbmerge` to obtain their
   marker tables.
3. With a [custom script](master/script/count_markers.pl), we then summarized the 
   [available markers](data/markers/markers.tsv). This indicated that the cpDNA markers 
   _rbcL_, _matK_, _rps16_ and _ndhF_ are suitable candidates.
4. For these candidate markers we then created a [table](data/markers/taxa.tsv) to collect 
   accession numbers for each of the species from step 2, using another 
   [custom script](script/count_taxa.pl).
5. Since the previous steps were based on an initial binning of the input taxa by 
   `taxa_in_original_tree` we may have filtered out putative PhyLoTA clusters that
   straddle multiple of these bins, which may have resulted in "missing" taxa. To correct
   this we then collected all the taxa for which we didn't yet have enough data and re-ran
   steps 1..4 on this taxon set. This resulted in a 
   [more enriched marker/taxa table](data/markers/merged.tsv).
   
## Additional enrichment
We performed more targeted searches for the identified markers to see if additional
sequence data are available (the SUPERSMART pipeline operates on an indexed version of
GenBank release 194, which is not the most recent one). This resulted in a
[final selection](data/enriched/final_selection_table_asterids.tsv).

## Multiple sequence alignment and phylogenetic inference
We then downloaded the sequences (using another [script](script/fetch_seqs.pl)) and aligned
each marker using MAFFT v7.130b (2013/12/05). We subsequently [merged](script/merge_aln.pl)
the alignments into a [supermatrix](data/enriched/merged.phy), which we analyzed using
ExaBayes 1.4.1 (by way of `smrt bbinfer`). We [rooted](script/reroot.pl) the trees in the 
posterior sample on the taxon _Cornales_.

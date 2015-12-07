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
   column of the [sampling_asterids_tree](data/taxa/sample_asterids_tree.tsv) spreadsheet,
   binned by their grouping under `taxa_in_original_tree`. For example, for the family
   Fouquieriaceae we expanded the genus Fouquieria.
2. For each of these `taxa_in_original_tree` we then ran `smrt align` and `smrt orthologize`
   to obtain their available PhyLoTA clusters, and then ran `smrt bbmerge` to obtain their
   marker tables.
3. With a custom script, we then summarized the [available markers](data/markers/markers.tsv).
   This indicated that the cpDNA markers _rbcL_, _matK_, _rps16_ and _ndhF_ are suitable
   candidates.
4. For these candidate markers we then created a [table](data/markers/taxa.tsv) to collect 
   accession numbers for each of the species from step 2.
   
## Additional enrichment
We then performed more targeted searches for the identified markers to see if additional
sequence data is available (the SUPERSMART pipeline operates on an indexed version of
GenBank release 194, which is not the most recent one). This resulted in a
[final selection](data/enriched/final_selection_table_asterids.tsv).

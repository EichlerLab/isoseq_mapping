# ISO-seq alignment pipeline.

Pipeline for generating ISO-seq alignments to a given reference genome.

## Adding a new collection of ISO-seq reads.

Add a new entry to the `config.json` file. It must contain the keyword `isoseq` to specify the input reads for that collection/species/group.

## Specify a reference genome.

Specify which reference genome you want to use by altering the `reference` value at the top level of the `config.json` file. It is important that any entry marked as reference *must* contain the following keywords: `reference`, `exons`, and optionally `gmap_db` and `gmap_db_name` and/or `star_index` depending on which alignment software you are using.

## Description of output files.

By default, the snakemake pipeline will generate a `{species}.gmap.filtered.sort.split.novel_exons.tbl` file. This is a BED file that contains novel exons (i.e., exons that do not occur in the reference transcripts) along with percent covered by segmental duplications and tandem repeats. BAM and BED files of the filtered alignments are stored in `{species}.gmap.filtered.bed`.

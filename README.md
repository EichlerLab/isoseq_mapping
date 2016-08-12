# ISO-seq alignment pipeline.

Pipeline for generating ISO-seq alignments to a given reference genome.

## Adding a new collection of ISO-seq reads.

Add a new entry to the `config.json` file. It must contain the keyword `isoseq` to specify the input reads for that collection/species/group.

## Specify a reference genome.

Specify which reference genome you want to use by altering the `reference` value at the top level of the `config.json` file. It is important that any entry marked as reference *must* contain the following keywords: `reference`, `exons`, and optionally `gmap_db` and `gmap_db_name` and/or `star_index` depending on which alignment software you are using.

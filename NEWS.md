# RiboseQC 1.2.0

## Installation and dependencies

- RiboseQC installs and runs again on current Bioconductor (3.22). Before,
  `prepare_annotation_files()` failed there, because GenomicFeatures'
  `makeTxDbFromGFF()` is defunct; it now comes from txdbmaker.
- devtools is no longer needed. `prepare_annotation_files(forge_BSgenome =
  TRUE)` installs the forged BSgenome package with `install.packages()`, and
  stops with an error if the installation fails.

## Changes in results

- `RiboseQC_analysis()` now gives the same P-site offsets in each run on the
  same data. Before, the offsets of some read lengths could change from one
  run to the next: `calc_cutoffs_from_profiles()` and `choose_readlengths()`
  call `kmeans()`, which picks random start centres, and RiboseQC set no
  seed. Each `kmeans()` call now uses a fixed seed (666) and R's default kinds
  of random numbers, through the withr package (which ggplot2 already
  installs). The random numbers of your R session don't change. So the
  offsets of these read lengths, and the P-sites, bedgraphs, `*_for_ORFquant`
  files and report plots that use them, can differ from those of an earlier
  run. The offsets of read lengths with a clear frame signal don't change.
- `prepare_annotation_files()` now has the code of ORFquant's copy of the
  function, which has the fixes below (ORFquant will use RiboseQC's). On
  GENCODE GTFs, the annotation is the same as before, apart from the changes
  under "Annotation files". On other GTFs, the biotypes and gene names can
  change:
  - They are read by attribute name, one row per transcript, from any line of
    the transcript (for genes, also from a line of the gene): biotypes from
    `gene_biotype` or `gene_type` and `transcript_biotype` or
    `transcript_type`, gene names from `gene_name`, `gene_symbol`, `gene`
    (NCBI) or `ref_gene_name` (StringTie). Before, the columns were named by
    their position, so a GTF with both `gene_type` and `gene_biotype` (or both
    `gene_name` and `gene_symbol`) gave mixed-up columns or an error.
    Biotypes that some transcripts lack are now `no_type`, not `NA`.
  - The transcript biotype `mRNA` (NCBI) is read as `protein_coding`. The
    bundled Arabidopsis GTF has it: its transcripts in `trann` and in
    `table_gene_tx_IDs` change from `mRNA` to `protein_coding`.
    `RiboseQC_analysis()` reads only the gene biotypes, so its results on
    this GTF don't change.
  - Without biotypes (for example in UCSC's GTFs), transcripts with CDS lines
    and their genes are `protein_coding`. Before, they were `no_type`, so all
    exonic regions outside the CDS were `ncRNAs` and none were `ncIsof`. For
    such GTFs, the reads counted in these regions and biotypes change.

## Annotation files

- New `*_Rannot` files have the format of ORFquant's. The new element
  `genome_package` holds the name of the forged BSgenome package, or `NULL`
  with `genome_seq`. The element `genome` now holds the genome sequence: the
  BSgenome object, or the `FaFile_Circ` as before. Before, it held the package
  name or the `FaFile_Circ`.
- `exons_bins` in `*_Rannot` keeps all columns of `exonicParts()` (`tx_id`,
  `tx_name`, `gene_id`, `exon_id`, `exon_name` and `exon_rank`). Before, it
  had `gene_id`, `tx_name` and an empty `exonic_part`. RiboseQC doesn't use
  `exons_bins`.
- `load_annotation()` reads the `*_Rannot` files of all versions of RiboseQC
  and ORFquant. Its new argument `envir` says where it puts `GTF_annotation`
  and `genome_seq`; the default is the environment of the caller, as before.
- The other files that `prepare_annotation_files()` writes don't change,
  apart from the biotypes in `table_gene_tx_IDs` (see "Changes in results").

## Bug fixes

- `load_annotation()` no longer fails with `first argument has length > 1`
  on annotations made with a forged BSgenome package (#21).
- `prepare_annotation_files()` no longer fails on GTFs without `transcript`
  lines, such as the bundled Arabidopsis GTF, or with lines without
  `transcript_id`, such as the `gene` lines of GENCODE and Ensembl GTFs
  (#18, #21).
- `prepare_annotation_files()` no longer fails with `<n> elements in value to
  replace <m> elements` when the GTF has genes with exons on both strands or
  on more than one chromosome, for example the genes in the pseudoautosomal
  regions of UCSC's RefSeq GTF. The annotation's `genes` leaves these genes
  out, as before, and their UTR, intron and non-coding exon regions get no
  gene id.
- `prepare_annotation_files()` no longer fails when the annotation has
  exactly one coding transcript.
- `prepare_annotation_files(forge_BSgenome = TRUE)` now keeps the circular
  chromosomes circular in the forged BSgenome package when the genome has
  more than one, for example `ChrM` and `ChrC`. Before, the seed file lost
  them.
- `prepare_annotation_files()` now checks its input before the long steps,
  and stops with a message that says what is wrong:
  - Files that don't exist (the GTF file, and the FASTA or twobit file) are
    listed together, before anything is written.
  - A GTF file that is empty, has no exon lines, has exon lines without
    `transcript_id` or `gene_id`, or has no CDS lines.
  - When txdbmaker cannot read the GTF file, the message says that its
    chromosome names must be those of the genome sequence, with examples.
  - Without `genome_seq` and without `twobit_file`, it asks for one.
- `forge_BSgenome` now defaults to `TRUE`, as its help said. With
  `genome_seq`, the forging is cancelled, with a message. Calls that worked
  before don't change: without `genome_seq`, `forge_BSgenome = FALSE` failed.
- `prepare_annotation_files()` now closes the TxDb database when it no longer
  needs it (lcalviell/ORFquant#3).
- `prepare_annotation_files()` gives the message "N genes were dropped
  because they have exons located on both strands..." once, not twice. It
  also gives a message with the number of transcripts that the TxDb has no
  exons for, for example because of trans-splicing.
- `prepare_annotation_files()` gives its progress lines with `message()`, so
  `suppressMessages()` hides them. Before, it wrote them to the standard
  output. It no longer writes the `.fai` index of the FASTA file again in
  each run: Rsamtools makes the index when it is missing.

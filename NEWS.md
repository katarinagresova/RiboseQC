# RiboseQC (development version)

## Input checks

- `RiboseQC_analysis()` now checks its input before it reads the first BAM
  file, and stops with a message that says what is wrong. The runs on correct
  input give the same results. The checks are:
  - Input files that don't exist (`annotation_file`, `bam_files`) are listed
    together. Before, the error for a missing annotation file was `cannot open
    the connection` (only a warning gave the file name), and a missing BAM
    file stopped the run only after the analysis of the BAM files before it.
    `annotation_file` must have 1 file, or 1 for each BAM file, as the help
    says. Before, fewer files stopped the run at the first BAM file without
    one, and more files than BAM files were not used.
  - The FASTA file of each annotation must exist. The annotation refers to it
    by its path, so a FASTA file that was moved or deleted after
    `prepare_annotation_files()` stopped the run near the end of the analysis
    of the first BAM file (in the P-sites profiles), with `file(s) do not
    exist`.
  - Each BAM file must have a chromosome of its annotation, and the
    chromosomes that the two have must have the same lengths. These are the
    checks of `prepare_for_ORFquant()`. Before, a BAM file with other
    chromosome names (`1` for `Chr1`), or aligned to the transcripts, stopped
    the run with `1 elements in value to replace 0 elements`, and a BAM file
    of another genome assembly with `incompatible seqlengths`. The error of
    `prepare_for_ORFquant()` for a BAM file without a chromosome of the
    annotation now names the two files, and says that the reads must be
    aligned to the genome.
  - The directories of the output files (`dest_names`, and `report_file`
    with `create_report = TRUE`) must be writable. A directory that doesn't
    exist is made, as before, but now before the first BAM file. Before, a
    directory that could not be written to stopped the run at the first file
    that it wrote, with `cannot open the connection`, and the directory of the
    report stopped it at the report, after the analysis.
- A BAM file without mapped reads (or with only secondary alignments and
  duplicates, which are not used) stops the run at the start of its
  analysis, with an error that says so. Before, the error was `subscript out
  of bounds`.
- `prepare_annotation_files()` warns when the lines of the GTF file do not
  agree: when the lines of a transcript give different values of an
  attribute (for example two biotypes), or the transcripts of a gene give
  different gene names or biotypes. Then the order of the lines chooses the
  value (see the help), and the warning names some of these transcripts and
  genes. The annotation does not change.
- `prepare_annotation_files()` stops at once for a GTF file without CDS
  lines also when the file has `start_codon` or `stop_codon` lines: `The GTF
  file ... has no CDS lines`. Before, txdbmaker made a CDS of 3 nt from each
  stop codon, so the annotation was made, and ORFquant stopped later with
  `subscript contains out-of-bounds indices`. The check reads the GTF file
  until it finds an exon line and a CDS line, before the TxDb is made.
- `RiboseQC_analysis()` warns, before the first BAM file, when a BAM file has
  reads on chromosomes that its annotation does not have. The warning gives
  the number of alignments on them, from the index of the BAM file. For a
  BAM file without an index, it names the chromosomes of the BAM header that
  the annotation does not have. These reads are used (for example in the
  P-site and coverage tracks), but they are in no region of the annotation.
  The results do not change.
- `prepare_annotation_files()` stops with an error that says what is wrong
  when a CDS of the GTF file ends after the end of its chromosome (one that
  is not circular) in the genome sequence, for example when the GTF file and
  the FASTA file are of different genome assemblies. The error names a
  transcript, its chromosome and the length of the chromosome. Before, the
  error was `Out of bounds ranges detected which are not on circular
  chromosomes`.

## Changes in results

- `RiboseQC_analysis()` on a BAM file that it reads in one chunk (with at
  most `chunk_size` alignments) gives the 5' ends of the reads
  (`read_stats$reads_pos1` in `<dest_name>_results_RiboseQC_all`) without
  names, as with more chunks. Before, each range and its score had its read
  length as its name, so the read lengths of the top 50 positions (`len_adj`
  of `sequence_analysis`, in `<dest_name>_results_RiboseQC` and
  `<dest_name>_results_RiboseQC_all`) were, for example, `29.29` instead of
  `29`, and their scores and percentages had names too. The positions, their
  widths and sequences, and all other values don't change. The report does
  not show `len_adj`.

## Bug fixes

- `RiboseQC_analysis()` no longer stops with `invalid 'file' argument` when
  `genome_seq` is an `FaFile`. The argument has no effect, as before: the
  genome sequence is the FASTA file of the annotation. A path was ignored
  without a message. Now a warning says so when `genome_seq` is another file
  than the FASTA file of the annotation. To use another FASTA file, make the
  annotation with it. The help of `genome_seq` said that it is used instead
  of a BSgenome package; it now says that it is not used.
- `RiboseQC_analysis()` no longer stops when the BAM file has reads on
  chromosomes that the annotation does not have (for example scaffolds), and
  the BAM file has more than `chunk_size` alignments, so that it is read in
  more than one chunk. The error was `The following seqlevels are to be
  dropped but are currently in use` (#14 reported `duplicated sequence
  names`). These reads now count as in a run with one chunk: they give
  P-sites and coverage on these chromosomes, and they count in the
  percentages of `P_sites_calcs`. The runs that worked before give the same
  results.
- `RiboseQC_analysis()` gives the same results for a BAM file with more than
  `chunk_size` alignments (5,000,000 by default), which it reads in more than
  one chunk, as for the same BAM file in one chunk. Before, the merge of the
  chunks joined adjacent positions with the same number of 5' ends or P-sites
  into one range, and the 5' and P-site profiles lost the reads of a range
  that crossed the end of an exon, or the start or the end of a CDS. So the
  profiles, the percentages of the "Top 50 mapping positions" (too high) and,
  when two offsets were near in the choice, the offset of a read length (and
  with it the P-site files and `for_ORFquant`) depended on `chunk_size`. The
  results of BAM files with more than `chunk_size` alignments can change. The
  results of the other BAM files do not change.
- `prepare_annotation_files()` no longer stops on a GTF file without
  introns, for example a GTF file of chrM only, or of transcripts with one
  exon each. `RiboseQC_analysis()` and `prepare_for_ORFquant()` no longer
  stop on an annotation without splice junctions, such as the annotation of
  such a GTF file. The error was `1 elements in value to replace 0 elements`.
  The results of the other annotations do not change.
- `prepare_annotation_files()` no longer stops with `'id_or_name2' must be a
  single string` or `argument is of length zero` when the chromosomes that
  are not circular have no genetic code or more than one: for example, a GTF
  file of chrM only with chrM circular, or a genome with chrM and
  `circ_chroms` without chrM. To find out if the CDS include their stop
  codons (`stop_in_gtf`), it now translates each CDS with the genetic code of
  its chromosome, and it uses the CDS on circular chromosomes if all CDS are
  on them. Before, it used one genetic code for all chromosomes that are not
  circular, also for those without CDS. The annotations that were made
  before do not change.

## Documentation

- The help of `prepare_annotation_files()` now says that the order of the
  lines in the GTF file does not change the results if the lines agree, as in
  GENCODE and Ensembl files (lcalviell/ORFquant#10). If the lines of a
  transcript, or the transcripts of a gene, give different values (for
  example two gene names), the order chooses the value: the first line in
  the file gives it.
- `?RiboseQC_analysis` now says that soft-clipped bases (`S` in the CIGAR)
  are not counted (#12): the read length is the aligned length, and the
  offset counts from the first aligned base at the 5' end. So a 31-nt read
  with 2 soft-clipped bases counts as a 29-nt read. Alignments with soft
  clipping (for example STAR's default, `--alignEndsType Local`) can give
  other read lengths, offsets and P-sites than end-to-end alignments.
  RiboseQC has always worked this way. `?prepare_for_ORFquant` and the
  README already said it.
- `?RiboseQC_analysis` now says which alignments are not used: alignments
  flagged as duplicates, secondary alignments and alignments with an
  insertion or a deletion (`I` or `D` in the CIGAR). These alignments are in
  no output: not in the read counts, the read length distributions, the
  offsets, the P-sites, the coverage or the junction reads, and the report
  does not count them. RiboseQC has always worked this way.
  `?prepare_for_ORFquant` and the README already said it. In the example
  data, 24 of the 26,702 alignments of the root sample and 26 of the 53,142
  alignments of the shoot sample have an insertion or a deletion. In a human
  sample aligned with STAR (SRR15513199), 0.7% of the alignments have one,
  but only 0.035% of the uniquely mapping alignments.
- The help of `prepare_annotation_files()` now says which annotations give
  wrong ORF categories in ORFquant. Since ORFquant 1.3.1, only an annotation
  with `stop_in_gtf` `NA` does (the codon after each CDS is its stop codon,
  also for a CDS that includes it). Before, the help said that a GTF file
  with both kinds of CDS gives wrong categories.

# RiboseQC 1.4.0

## Installation

- The new file `environment.yml` describes a conda environment with R 4.5,
  Bioconductor 3.22, the packages that RiboseQC needs and pandoc. The README
  shows how to install RiboseQC into it from GitHub, or with BiocManager.

## New arguments

- `RiboseQC_analysis()` has the new argument `n_cores` (default 1). With
  `n_cores` more than 1, it reads the BAM file in chunks of `chunk_size`
  alignments, as before, and analyses the chunks in forked processes, at most
  `n_cores` at a time. Each process needs the memory for one chunk. The
  results are the same as with 1 core. Forked processes are not available
  on Windows.

## Changes in results

- `RiboseQC_analysis()` and `prepare_for_ORFquant()` find the uniquely
  mapping reads with mismatches (`P_sites_uniq_mm` in
  `<dest_name>_for_ORFquant`, and in `<dest_name>_results_RiboseQC_all`) by
  the mismatches in the MD tag: a number and then a base. Before, a read had
  mismatches when its MD tag had more than 3 characters. So a read shorter
  than 20 nt with 1 mismatch (for example `5A9`) did not count, and a read of
  1,000 nt or more without mismatches did. Only these reads change. A read
  shorter than 20 nt gets a P-site only when its read length has a cutoff:
  in `RiboseQC_analysis()`, for example with `rescue_all_rls = TRUE` or
  `offsets_df`, and in `prepare_for_ORFquant()`, in the cutoff file. BAM
  files without MD tags give no `P_sites_uniq_mm`, as before. In ORFquant,
  these P-sites give only the column `P_sites_raw_uniq_mm`. The other P-sites
  and the other output files don't change.

## Bug fixes

- `prepare_annotation_files()` stores the path of the FASTA file
  (`genome_seq`) as an absolute path. Before, it stored the path as given, so
  with a relative path, `load_annotation()`, `RiboseQC_analysis()` and
  ORFquant failed in another working directory: the FASTA file did not
  exist. Only the directories of the path are resolved, not the file name, so
  the `.fai` and `.gzi` files next to the FASTA file are used, also when the
  file is a symbolic link. The path of an `FaFile` is stored as it is.

## Documentation

- The new vignette, `vignette("RiboseQC")`, runs the analysis on the example
  data included in the package, and shows some of the results: the reads of
  each compartment, the read lengths, the P-site offsets and the frame of the
  P-sites. It also shows that `prepare_for_ORFquant()` gives the same P-sites
  as `RiboseQC_analysis()`. To install it, see the README.
- The README now installs RiboseQC from this repository (#7, #19). Before,
  it installed lcalviell/Ribo-seQC, which has older code. It also runs the
  analysis on the example data included in the package, and describes the
  input and output files.
- The repository no longer holds `RiboseQC-manual.pdf`. It was the manual of
  version 0.99.0 from February 2019, without the changes to the help pages
  since then. Read the help pages with `help(package = "RiboseQC")`.
- The old vignette, `RiboseQC.Rmd`, and its code, `RiboseQC.R`, are
  removed. The vignette was not installed with the package. It downloaded
  its data from a server that no longer answers, and forged BSgenome
  packages, which RiboseQC no longer does.
- The example data no longer include `shoots_results_RiboseQC`: no example
  used it. The repository no longer holds `ex/`, the output files of an old
  analysis of the example data. The removed files are still in the
  repository's history, for example at commit `7d2eeac`.
- The new file `README.md` of the example data (in
  `system.file("ext_data", package = "RiboseQC")`) gives the sources of the
  files, their terms of use and the papers to cite.
- The genome sequence of the example data, `simp_arab.fa.gz`, is compressed
  with `bgzip` (before: gzip) and comes with its indexes (`.fai`, `.gzi`). The
  sequence is the same. `prepare_annotation_files()` reads it directly, so the
  README, the vignette and the examples no longer write an uncompressed copy
  first. Code that reads the file with gzip tools, for example
  `Biostrings::readDNAStringSet()`, still works.
- The example file `root_results_RiboseQC` is new: the results of
  `RiboseQC_analysis()` with its default parameters on `simp_arab_root.bam`.
  The old file came from an earlier and larger subset of the root sample, and
  from the code of 2019. The examples of the help pages use this file, so
  their plots change. `inst/script/make_root_results.R` makes the file.
- The example of `plot_codon_usage_bulk_rmd()` runs: on the new
  `root_results_RiboseQC`, it takes about 12 s (before: about 1.5 min, in
  `\dontrun{}`).
- Corrections to the help pages:
  - `?RiboseQC_analysis`: `chunk_size` must be more than 100,000 (before:
    10,000). The output files get the prefixes of `dest_names`, which default
    to the BAM paths (before: the BAM paths). The report is by default the
    first BAM path followed by `_RiboseQC_report.html` (before: `.html`).
    `offsets_df` can be a data.frame, for the offsets of `nucl`, or a list with
    a table for each compartment.
  - `?create_html_report`: the plots are saved as RDS files in
    `<output_file>_plots/rds/` (before: as PDF files, in `<name>_plots/`).
    `create_pdfs_from_rds_objects()` makes PDF files of them.

# RiboseQC 1.3.0

## New functions

- `prepare_for_ORFquant()` comes from ORFquant, with the same code and the
  same results. It calculates the P-sites and the junction reads of a BAM
  file with a cutoff table, for example `<dest_name>_P_sites_offsets.tsv` of
  `RiboseQC_analysis()`, and writes the file `<dest_name>_for_ORFquant` that
  ORFquant's `run_ORFquant()` reads. With `n_cores` more than 1, it reads the
  BAM file in chunks and calculates the P-sites of the chunks in forked
  processes. ORFquant imports it from RiboseQC and exports it again, so
  scripts that call it after `library(ORFquant)` don't change.
- RiboseQC exports `check_n_cores()`, `check_output_dir()`, `mc_job()`,
  `mc_value()` and `mc_lapply()` for ORFquant, which uses them too. They are
  not for users, and their help pages are not in the index.

# RiboseQC 1.2.1

## New output file

- `RiboseQC_analysis()` writes the new file `<dest_name>_P_sites_offsets.tsv`.
  It has the offsets that it used to calculate the P-sites, with the columns
  `read_length`, `cutoff` and `compartment` (`nucl` or the name of a circular
  chromosome). It includes the offsets of `offsets_df`, and with
  `rescue_all_rls = TRUE` the offset 12 of the other read lengths. It is the
  cutoff file of ORFquant's `prepare_for_ORFquant()`. Use it, not
  `*_P_sites_calcs`: that file lists the offsets of the analysis for all
  read lengths, and not only those that were used for the P-sites.

## Changes in results

- `RiboseQC_analysis()` with `normalize_cov = TRUE` (the default) no longer
  stops at "Exporting files" with `Scores must be non-NA numeric values`
  (#15). This error came when the coverage summed to less than 1 million,
  for example on the bundled root BAM. The normalised coverage tracks
  (`*_coverage_plus.bedgraph`, `*_coverage_minus.bedgraph`,
  `*_coverage_uniq_plus.bedgraph` and `*_coverage_uniq_minus.bedgraph`) now
  sum to 1 million, as documented. Before, the function divided the integer
  coverage by the correction factor (the sum of the coverage / 1 million),
  and IRanges first rounded this factor down to an integer. A factor below 1
  gave a division by 0. A factor of 1 or more gave values that were too
  large, by the factor divided by its integer part (less than 2 times). So
  the values in these four files change. The tracks of
  `normalize_cov = FALSE` and the other output files don't change.
- `RiboseQC_analysis()` calculates the P-sites of uniquely mapping reads with
  mismatches (`P_sites_uniq_mm` in `<dest_name>_for_ORFquant` and in
  `<dest_name>_results_RiboseQC_all`) as ORFquant's `prepare_for_ORFquant()`
  does. With the offsets of `<dest_name>_P_sites_offsets.tsv`, the two
  functions now give the same `P_sites_uniq_mm`. Before, the values were
  wrong. On the minus strand, the function kept the reads and not their
  P-sites, so each position of a read counted. On the plus strand, it
  selected the P-sites of other reads. A read has mismatches when its MD tag
  has more than 3 characters. Before, 3 characters were enough (for example
  `5A5`). This difference applies only to reads shorter than 20 nt and reads
  of 100 nt or more. BAM files without MD tags give no `P_sites_uniq_mm`, as
  before. In ORFquant, these P-sites give only the column
  `P_sites_raw_uniq_mm`. The other P-sites and the other output files don't
  change.

## Bug fixes

- The HTML report (`create_html_report()`, and `RiboseQC_analysis()` with
  `create_report = TRUE`, the default) works again with current knitr
  (tested with 1.52). Before, rendering stopped with an error about a code
  chunk, because the report templates indented the closing fences of the
  chunks that they generate.
- After a report fails, R output goes to the console again. Before,
  `create_html_report()` left its text-output file open, so R wrote all
  later output of the session into it.
- When `RiboseQC_analysis()` skips the P-sites and the codon occupancy (not
  enough signal or a low frame preference), it now gives a warning, not a
  printed message. The message now names the correct argument,
  `readlength_choice_method = 'all'`. Before, it named `choose_readlengths`.

# RiboseQC 1.2.0

## Installation and dependencies

- RiboseQC installs and runs again on current Bioconductor (3.22). Before,
  `prepare_annotation_files()` failed there, because GenomicFeatures'
  `makeTxDbFromGFF()` is defunct; it now comes from txdbmaker.
- RiboseQC now needs Bioconductor 3.22 (R 4.5) or later. It uses the
  cigarillo package, new in 3.22, in place of GenomicAlignments'
  `cigarRangesAlongQuerySpace()` and `cigarRangesAlongReferenceSpace()`,
  which are deprecated. So `RiboseQC_analysis()`, `get_ps_fromspliceplus()`
  and `get_ps_fromsplicemin()` no longer give deprecation warnings.
- devtools and GenomicFiles are no longer needed.
- `prepare_annotation_files()` no longer forges and installs a BSgenome
  package from a 2bit file: this failed with current BSgenomeForge versions
  (#21, #23). Give the genome as a FASTA file (`genome_seq`); `twoBitToFa` of
  the UCSC tools makes one from a 2bit file. `forge_BSgenome` defaults to
  `FALSE`, as before, and without `genome_seq` the function stops at once
  with this advice. `twobit_file`, `scientific_name` and `annotation_name`
  are no longer used. Annotations that earlier versions made with a forged
  BSgenome package still load.
- `library(RiboseQC)` now attaches only GenomicRanges and the packages it
  attaches (IRanges, S4Vectors, Seqinfo, BiocGenerics, generics and stats4).
  Before, it attached the 18 packages of its Depends field and theirs. They
  are now imported, so RiboseQC still uses them. Scripts that relied on
  `library(RiboseQC)` to attach others, for example Biostrings,
  rtracklayer, Rsamtools or ggplot2, need their own `library()` calls. The
  vignette now calls `Biostrings::getGeneticCode()` and
  `ggpubr::as_ggplot()`. jsonlite is no longer needed: RiboseQC didn't use
  it.
- Loading RiboseQC no longer warns "replacing previous import" for
  `combine` and `Position`.
- A new run of roxygen2 on the source no longer breaks the install (#24).
  Before, it added an import of `DEFAULT_CIRC_SEQS` from GenomicFeatures,
  which no longer has it. This list of circular chromosome names (the
  default of `circ_chroms` in `prepare_annotation_files()`) is now a
  constant of RiboseQC, with the same names. `NAMESPACE` and the help pages
  are now made with roxygen2 from the source.

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
- `prepare_annotation_files()` now chooses the representative transcript of
  each gene from that gene's own start and stop codons. Before, when the
  start and stop codons of two genes alternated along a chromosome (for
  example in overlapping genes), the flags for a gene's most common and
  outermost codons could go to the other gene's codons. So `upstr_stasto`,
  `mostfreq_stasto` and the columns `reprentative_mostcommon`,
  `reprentative_boundaries` and `reprentative_5len` were wrong for some
  transcripts and genes: on GENCODE v47, 339 of 20,503 genes get another
  `reprentative_mostcommon` transcript. `RiboseQC_analysis()` takes these
  transcripts for the metagene of the nuclear CDS, so with a new annotation
  its 5' profiles, P-site offsets, P-sites, bedgraphs, `*_for_ORFquant` files,
  results files and report plots can change. The annotation of the bundled
  example doesn't change. Annotations made by earlier versions keep their
  values, and give the same results as before.

## Annotation files

- New `*_Rannot` files have the format of ORFquant's: they have the new
  element `genome_package` (`NULL`), and the element `genome`, which holds
  the `FaFile_Circ` as before, is now the last one.
- `exons_bins` in `*_Rannot` keeps all columns of `exonicParts()` (`tx_id`,
  `tx_name`, `gene_id`, `exon_id`, `exon_name` and `exon_rank`). Before, it
  had `gene_id`, `tx_name` and an empty `exonic_part`. RiboseQC doesn't use
  `exons_bins`.
- `load_annotation()` reads the `*_Rannot` files of all versions of RiboseQC
  and ORFquant. Its new argument `envir` says where it puts `GTF_annotation`
  and `genome_seq`; the default is the environment of the caller, as before.
- `load_annotation()` now returns the annotation and the genome sequence,
  invisibly, as a list with the elements `GTF_annotation` and `genome_seq`,
  as ORFquant's `load_annotation()` does. It still assigns them in `envir`,
  so scripts that use these variables work unchanged. Before, it returned
  `GTF_annotation` alone, invisibly, which was not documented: a script that
  used this value must now take `$GTF_annotation` of the list.
- The other files that `prepare_annotation_files()` writes don't change,
  apart from the biotypes in `table_gene_tx_IDs` and the columns
  `upstr_stasto`, `mostfreq_stasto` and `reprentative_*` in
  `cds_txs_coords_similbed.bed` (see "Changes in results").

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
- `prepare_annotation_files()` now checks its input before the long steps,
  and stops with a message that says what is wrong:
  - Files that don't exist (the GTF and FASTA files) are listed together,
    before anything is written.
  - A GTF file that is empty, has no exon lines, has exon lines without
    `transcript_id` or `gene_id`, or has no CDS lines.
  - When txdbmaker cannot read the GTF file, the message says that its
    chromosome names must be those of the genome sequence, with examples.
  - Without `genome_seq`, it asks for a FASTA file (see "Installation and
    dependencies"). Before, it failed in `FaFile(NULL)`.
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
- The `getSeq()` error for ranges that wrap twice around a circular
  chromosome now reads "Ranges wrapping twice isn't implemented yet...",
  without the stray quote, line break and spaces it had.
- RiboseQC now imports all the functions that it uses from other packages,
  for example `Seqinfo()`, `isCircular()` and `saveDb()`. Before, R found
  some of them only because `library(RiboseQC)` attached the packages in
  `Depends`. When RiboseQC was only loaded (with `RiboseQC::`, with
  `loadNamespace()` or as an import of another package), functions such as
  `seqinfo()` of a `FaFile_Circ` stopped with `could not find function
  "Seqinfo"`.
- The help pages are made again from the source (with roxygen2 8.1.0; they
  were from 2019). They now describe all arguments, for example
  `genome_seq`, `stranded`, `normalize_cov` and `offsets_df` of
  `RiboseQC_analysis()`, and their usage sections match the functions. The
  examples now run on the bundled data. Before, they used data that the
  package doesn't have (`data(res_all)`) or files that don't exist. Four of
  them are not run: those of `plot_codon_usage_positional_rmd()` and
  `plot_codon_usage_bulk_rmd()` take minutes, and those of
  `calc_cutoffs_from_profiles()` and `choose_readlengths()` need a
  `*_results_RiboseQC_all` file.
- `?generate_rdata_list` is now installed with the package. Before,
  `R CMD build`, which `remotes::install_github()` also runs, left it out.
- `?prepare_annotation_files` now says what `stop_in_gtf` holds (`"*"` if
  most CDS end with a stop codon, `NA` if most end before it) and how
  `stop_codon` lines are used. The help of `getSeq()` for `FaFile_Circ` now
  says how it reads ranges past the end of a circular chromosome, and what it
  returns.

# RiboseQC 0.99.0

- Submitted to Bioconductor (2018-05-15).

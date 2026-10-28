# RiboseQC (development version)

## Installation and dependencies

- RiboseQC no longer has methods of `unlist()` for `GRanges` and
  `GAlignments`. These methods made `unlist(x)` give `x` back, in the whole R
  session, when RiboseQC or ORFquant was loaded. So they also changed
  `unlist()` in other packages and in scripts: without them, `unlist()` of a
  `GRanges` stops with an error, and `unlist()` of an empty `GRanges` gives
  `NULL`. The results of RiboseQC do not change. `run_ORFquant()` of ORFquant
  1.3.1 and older needs the method for `GRanges`, so it stops with an error
  with this version of RiboseQC. Use ORFquant 1.3.2 or newer.

## New features

- The script `run_riboseqc.R`, installed with the package
  (`system.file("scripts", "run_riboseqc.R", package = "RiboseQC")`), runs
  `prepare_annotation_files()` and `RiboseQC_analysis()` from the command
  line, each in a new R process, for example `Rscript run_riboseqc.R --gtf
  genes.gtf --fasta genome.fa --bam root.bam,shoots.bam --outdir results`,
  and lists the files it writes. `--annotation` reuses the annotation of an
  earlier run; `--sample-names`, `--stranded`, `--readlength-method` and
  `--cores` are passed to `RiboseQC_analysis()`; `--no-report` and
  `--extended-report` choose the report; the other parameters have their
  default values. It checks its options before the first step. `--help`
  lists all options.
- The script `run_riboseqc.sbatch`, installed next to `run_riboseqc.R`, runs
  it as a SLURM job: `sbatch run_riboseqc.sbatch` followed by
  `run_riboseqc.R`'s options, with `--cores` set to the job's CPUs. It asks
  for 1 CPU, 32 GB of memory and 12 hours, and its comments give the time
  and memory measured on human samples.
- The progress messages of `RiboseQC_analysis()`, such as "Analyzing BAM
  file: ... ", now come from `message()`, not from `cat()`, as the messages
  of `prepare_annotation_files()` and `prepare_for_ORFquant()`. So they go
  to the standard error, not to the standard output, and
  `suppressMessages()` hides them. Their text and date are the same. The
  scripts `run_riboseqc.R` and `run_riboseqc.sbatch` show them as before.

## New arguments

- `RiboseQC_analysis()` has the new argument `write_TSV_files`, as
  `run_ORFquant(write_TSV_file = )` of ORFquant. With `TRUE`, the default,
  it writes the main tables of the results of each BAM file as
  tab-separated files with the prefix `<dest_name>`:
  `_read_lengths.tsv` (the reads of each read length in each compartment),
  `_regions.tsv` (the reads on each region, for example the CDS, in each
  compartment), `_genes.tsv` (the reads, RPKM and TPM of each gene),
  `_P_sites_calcs.tsv` (the same table as `_P_sites_calcs`),
  `_top50_positions.tsv` (the 50 positions with the most reads, with their
  sequences) and `_summary.tsv` (for each compartment: the reads, the
  percentage of them on the CDS, the percentage of them that give P-sites,
  and the frame preference of these reads). `?RiboseQC_output` gives
  their columns. The other output files do not change. With `FALSE`, the
  function writes none of these files, as before. A value other than `TRUE`
  or `FALSE` stops the function before it reads a file.
- `prepare_annotation_files()` has the new argument `genetic_codes`: the
  genetic codes of some chromosomes, e.g. `genetic_codes = c(chrM = "3")`,
  with the ids of `Biostrings::GENETIC_CODE_TABLE`. They replace the codes
  that the chromosome names give: `chrM`, `MT`, `MtDNA`, `mit` and
  `mitochondrion` get code 2 (vertebrate mitochondria), `Mito` gets code 3
  (yeast mitochondria), `dmel_mitochondrion_genome` gets code 5 (invertebrate
  mitochondria), and all other chromosomes get code 1. So the mitochondrion
  of a yeast or Drosophila genome with UCSC names (`chrM`) gets code 2. Then
  the codon tables of `RiboseQC_analysis()` give the amino acids of code 2
  for its codons, and ORFquant finds and translates its ORFs with code 2.
  The default, `NULL`, gives the same codes as before. A code that is not an
  id, or a chromosome that is not in the genome, stops the function before
  it makes the TxDb. The help of `prepare_annotation_files()` now gives these
  rules.
- `prepare_for_ORFquant()` has the new argument `stranded`, with the values
  of `RiboseQC_analysis()`: `TRUE` (the default), `FALSE` or `"inverse"`. Use
  `"inverse"` for libraries with the reads on the strand opposite to the RNA:
  the function then inverts the strand of each read before it calculates the
  P-sites and the junction reads, as `RiboseQC_analysis(stranded =
  "inverse")` does. Before, it used the strands of the BAM file, so for these
  libraries each P-site was on the other strand, with its offset counted from
  the 3' end of the read, and the reads on the junctions were not counted.
  `TRUE` and `FALSE` give the same results as before. Another value stops the
  function before it reads a file.

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
- `load_annotation()`, and so `RiboseQC_analysis()` and
  `prepare_for_ORFquant()`, stop with an error that says what is wrong when
  the annotation file is not a file of `prepare_annotation_files()`, for
  example the GTF file. Before, the error was `bad restore file magic
  number`.
- `RiboseQC_analysis()` stops at once when a value of `stranded` is not
  `TRUE`, `FALSE` or `"inverse"`: `'stranded' must be TRUE, FALSE or
  'inverse', not reverse`. Before, a value such as `"reverse"` or `"yes"`
  stopped the run in the analysis of the first BAM file, with `'ignore.strand'
  must be TRUE or FALSE`, and `NA` with `missing value where TRUE/FALSE
  needed`; a large BAM file was read once before. The values that worked,
  also `"TRUE"` or `1`, still work.
- `prepare_for_ORFquant()` stops at once when a compartment of the cutoff
  file is a chromosome that is not circular, for example `chr22`: a
  compartment is `nucl` (all chromosomes that are not circular) or the name
  of a circular chromosome. Before, the reads of that chromosome got P-sites
  with the rows of the chromosome and also with the rows of `nucl`, so a
  read length in the two gave each of its reads two P-sites. The help of
  `path_to_rl_cutoff_file` now says what a compartment is.

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
- `RiboseQC_analysis()` on a BAM file that it reads in more than one chunk
  (with more than `chunk_size` alignments) gives the coverage of the reads
  on the minus strand (`P_sites_stats` in
  `<dest_name>_results_RiboseQC_all`) the names `coverage_all_minus` and
  `coverage_uniq_minus`, as with one chunk and as `?RiboseQC_output` says.
  Before, their names were `coverage_all_min` and `coverage_uniq_min`. Their
  values, and all other output files, don't change. Code that reads them
  with `$coverage_all_min` still gets them, because `$` also finds a name by
  its start. Code that uses `[["coverage_all_min"]]` must use the new names.
- `RiboseQC_analysis()` on an annotation with fewer than 250 transcripts:
  the P-sites profiles and the codon counts (`profiles_P_sites`) now use the
  transcripts that the 5' profiles use (one representative transcript of
  each gene, with a CDS of more than 100 nt). Before, when 5 or more of
  these transcripts had P-sites, they used all transcripts of the
  annotation, so a transcript without a CDS stopped the run with `subscript
  contains invalid names`, and a transcript with a short UTR or CDS with
  `only 0's may be mixed with negative subscripts`: for example the
  annotation of a GTF file of chrM. The results of larger annotations do
  not change.
- `RiboseQC_analysis()` with circular chromosomes: the `nucl` row of
  `read_stats$rld_unq` (in `<dest_name>_results_RiboseQC` and
  `<dest_name>_results_RiboseQC_all`) now counts only the uniquely mapping
  reads of the other chromosomes, as the `nucl` row of `read_stats$rld` does
  for all reads. Before, it also counted the uniquely mapping reads of the
  circular chromosomes, which also have their own rows. For example, the
  `nucl` row of the bundled root sample had 22,095 reads, and now has
  17,154: without the 1,768 reads of ChrC and the 3,173 reads of ChrM. All
  other values, and the results of annotations without circular
  chromosomes, do not change. The report does not use `rld_unq`.
- `RiboseQC_analysis()` and `prepare_for_ORFquant()`: the P-sites of the
  uniquely mapping reads (`P_sites_uniq`) and of the uniquely mapping reads
  with mismatches (`P_sites_uniq_mm`) now come from the correct reads.
  Before, a read with 2 or more introns and its P-site in a middle block
  could get the mapping quality and the MD tag of another such read of the
  same read length and strand. So these P-sites could include the P-site of
  a read with a mapping quality of 50 or less, and not the P-site of a read
  with a higher mapping quality. This changes `P_sites_uniq` and
  `P_sites_uniq_mm` in `<dest_name>_for_ORFquant` and in `P_sites_stats`
  of `<dest_name>_results_RiboseQC_all`, and the files
  `<dest_name>_P_sites_uniq_plus.bedgraph` and
  `<dest_name>_P_sites_uniq_minus.bedgraph`. The P-sites of all reads
  (`P_sites_all`) do not change. For this fix, `get_ps_fromspliceplus()`
  and `get_ps_fromsplicemin()` now give the P-sites in the order of the
  reads, not sorted by position.
- `RiboseQC_analysis(stranded = FALSE)`: the reads of each region by read
  length (`read_stats$reads_summary` and `reads_summary_unq`, in
  `<dest_name>_results_RiboseQC` and `<dest_name>_results_RiboseQC_all`) now
  ignore the strand of the regions, as the reads of each region
  (`read_stats$positions`) and of each gene (`counts_*`) do. So the sum of
  each row is the value of `positions`. Before, a read counted only for the
  regions on its strand. For example, a read on the CDS of a gene on the
  plus strand now also counts for an intron of a gene on the minus strand
  at the same position. The report shows these counts in sections 1.2 and
  3, so they now agree with section 1.1 (`positions`). `read_subset`
  selects the read lengths from them, so with `stranded = FALSE` the run
  can now calculate the offsets of other read lengths, and its P-sites can
  change. The results of `stranded = TRUE` and `"inverse"` do not change.
- `RiboseQC_analysis()`: the 5' profile of `"all"` of a compartment
  (`profiles_fivepr$five_prime_bins` and `$five_prime_subcodon`, in
  `<dest_name>_results_RiboseQC_all`) is now the sum of the profiles of the
  read lengths, as `?RiboseQC_output` says, also when a transcript of the
  compartment has a UTR shorter than 50 nt. Before, then it was the profile
  of the last (longest) read length. For example, in the bundled root
  sample, `"all"` of ChrC was the profile of the 45 nt reads, and `"all"` of
  ChrM the profile of the 35 nt reads. `nucl` is affected when the
  transcripts used for it have such UTRs, for example with an annotation
  without UTRs. This changes the `"all"` row of
  `selection_cutoffs$results_cutoffs`, the `"all"` element of
  `selection_cutoffs$analysis_frame_cutoff`, and the plots of `"all"` in
  sections 4.1 and 4.2.1 of the report. For example, ChrC of the bundled
  root sample now has an `"all"` row (cutoff 21, frame preference 39), and
  before, it had none, because the profile of the 45 nt reads gave no
  offset. The profiles of the single read lengths, the offsets, the
  selected read lengths and the P-sites do not change.

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
- `RiboseQC_analysis(stranded = "inverse")` now gives the results of
  `stranded = TRUE` on the same reads with inverted strands. Before, only the
  offsets came from the inverted reads. The P-sites, the coverage tracks and
  the junctions came from the strands of the BAM file, so each P-site was on
  the other strand, with its offset counted from the 3' end of the read. This
  changed the P-site files, `for_ORFquant` and the P-site profiles and codon
  statistics of the report. The reads of each region and gene were counted
  on both strands, as with `stranded = FALSE`. The results of
  `stranded = TRUE` and `FALSE` do not change. For `prepare_for_ORFquant()`,
  see its new argument `stranded`.
- `getSeq()` of a `FaFile_Circ` now reads a range that starts after the end
  of a circular chromosome: the range is first moved back by a multiple of
  the chromosome's length, so that it starts on the chromosome. Before, such
  a range was dropped, so the sequences did not match the ranges (with the
  ranges `chrM:12-14` and `chr2:1-3` on a chrM of 10 nt, the result was one
  sequence, the sequence of `chr2:1-3` followed by that of `chrM:2-4`), or
  the function stopped with `subscript is a logical vector with out-of-bounds
  TRUE values`, or R crashed (a range that starts just after the end). The
  other ranges give the same sequences as before.
- The messages of `load_annotation()` now say where it puts
  `GTF_annotation` and `genome_seq`: "the global workspace" when `envir` is
  the global environment, else "the environment envir". Before, they said
  "the global workspace" and "parent workspace" for any `envir`, also in
  `RiboseQC_analysis()`, which loads the annotation into its own
  environment. The errors of `get_metagene_data()` and
  `get_codon_usage_data()` for a wrong `comp` now say `Variable comp`, not
  `Variable res`.
- `get_top50_cds_genes()` and `get_top50_all_genes()` give all genes when the
  results have fewer than 50 (100) genes. Before, they added rows of `NA` up
  to 50 (100) rows, which the "Top 50 abundant genes" tables of the report
  showed. The results with more genes give the same tables.

## Performance

- `RiboseQC_analysis()` is faster when it makes the P-site profiles for many
  read lengths, for example with `rescue_all_rls = TRUE` or
  `readlength_choice_method = "all"`. For the codon tables, it extracted the
  sequences of the same transcripts from the FASTA file once for each read
  length. It now extracts them once for each compartment. On a human sample
  (a 69 MB BAM, GENCODE 47, `rescue_all_rls = TRUE`, 32 read lengths), a run
  with the report on one core takes 392 s instead of 444 s (12% less). The
  results don't change.
- `RiboseQC_analysis()` is faster for BAM files that it reads in two or more
  chunks (more than `chunk_size` alignments, 5 million by default). It merged
  the 5' ends and the P-sites of each chunk with the sums of the chunks
  before it, so each merge took longer than the one before. It now sums the
  chunks once, after the last chunk. On a human sample (a 2.3 GB BAM, 18
  chunks, GENCODE 47, `rescue_all_rls = TRUE`), a run with the report on one
  core takes 2,403 s instead of 3,025 s (21% less), with the same memory. The
  results don't change.
- `RiboseQC_analysis()` and `prepare_for_ORFquant()` find the P-sites of each
  chunk of the BAM file faster. On a human sample (a 2.3 GB BAM, 18 chunks,
  GENCODE 47, `rescue_all_rls = TRUE`), a chunk of 5 million alignments takes
  27 s instead of 44 s, and a run of `RiboseQC_analysis()` with the report on
  one core takes 2,064 s instead of 2,393 s (14% less), with the same memory.
  The results don't change.
- The HTML report is faster: it made each heatmap with the scale "zscore"
  (the metagene heatmaps and the positional codon usage) two times. On a
  human sample (a 69 MB BAM, GENCODE 47, `rescue_all_rls = TRUE`, 32 read
  lengths), the report takes 136 s instead of 157 s. The plots don't change.
- `RiboseQC_analysis()` made the P-site profiles of all read lengths together
  ("all") two times for each compartment. It now makes them once. On the same
  sample (a 69 MB BAM), the P-site profiles take 90 s instead of 98 s. The
  results don't change.
- `RiboseQC_analysis()` reads the BAM file and counts the reads faster. It
  counted the reads of each read length in the regions (CDS, UTRs, introns,
  ...) with one `summarizeOverlaps()` for each read length, and it found the
  overlaps of the uniquely mapping reads a second time. It now finds the
  overlaps of all reads once. It also reads only the fields of the BAM file
  that each pass uses, and it finds the alignments with an insertion or a
  deletion and the read length without the soft-clipped bases faster, also in
  `prepare_for_ORFquant()`. On a human sample (a 2.3 GB BAM, 18 chunks,
  GENCODE 47, `rescue_all_rls = TRUE`), the counts take 327 s instead of
  506 s, the read lengths 95 s instead of 139 s, and a run with the report on
  one core takes 1,750 s instead of 2,000 s (12% less), without more memory.
  The results don't change.

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
- Corrections to the help pages. The code does not change.
  - `?RiboseQC_analysis`: the value of the function is now in the help: the
    paths of the `_results_RiboseQC` files, one for each BAM file. The help
    lists the output files and when the function writes them:
    `_results_RiboseQC_all` only with `write_tmp_files = TRUE`,
    `_for_ORFquant` and the P-site bedgraph files only if there are P-sites,
    the PDF files of `pdf_plots` only with the report. `write_tmp_files` is
    one value for all BAM files (before: also one for each BAM file, which
    stops the run). `read_subset` selects the read lengths by their reads on
    the CDS and the UTRs. `rescue_all_rls` gives the offset 12 to all read
    lengths without an offset (before: those without enough reads).
    `fast_mode` uses 500 transcripts (before: genes). `offsets_df` changes
    the P-sites and what the function makes from them, not the offsets in
    `selection_cutoffs`, `summary_P_sites`, `_P_sites_calcs` and the report.
    In the description of the results: `counts_cds_genes`,
    `counts_all_genes` and `reads_summary` are in `read_stats`; the P-sites
    of the uniquely mapping reads are `P_sites_uniq` (before:
    `P_sites_all_uniq`), and these reads have a mapping quality of more than
    50; `sequence_analysis` is a GRanges; `Codon_counts` is summed over the
    transcripts; the first window of 33 nt starts with the start codon
    (before: after it), and the last one ends at the end of the CDS (before:
    before the stop codon). The help now names the `_unq` statistics, the A-site
    and E-site tables and the smaller `_results_RiboseQC` file.
  - `?prepare_annotation_files`: the value of the function is the path of
    the `_Rannot` file (before: the help said only that it makes files), and
    the help lists the files that it makes. `annot_file` is made with
    `save()` (before: an "rds file"). The help lists the 20 default names of
    `circ_chroms` (before: `DEFAULT_CIRC_SEQS`, which is not exported). The
    CDS of transcripts with a CDS shorter than 3 nt are not used, but the
    transcripts stay in the annotation (before: "are removed"). The columns
    of `cds_txs_coords` and `trann` are as in the annotation.
  - `?load_annotation`: only the annotations of ORFquant with
    `forge_BSgenome = TRUE` have a BSgenome package in `genome_package`
    (before: also those of RiboseQC 1.2.0 and later, which never have one).
  - `?getSeq` of `FaFile_Circ`: a range must start on the circular
    chromosome to run past its end.
  - `?prepare_for_ORFquant`: the value of the function is now in the help:
    invisibly, the path of the `_for_ORFquant` file. The help says that the
    function assigns `GTF_annotation` and `genome_seq` in the global
    environment, and which column names the cutoff file needs.
    `?forked_processes`, `?input_checks`: `RiboseQC_analysis()` uses these
    functions too (before: only `prepare_for_ORFquant()`).
  - The help pages of the report functions. The plot functions print their
    plot and return `NULL` (before: "returns a plot"); the codon usage plots
    return the plot. The help gives the name of the RDS file that each one
    saves, and `output_rds_path` must end with `/`. The `_rmd` functions save
    RDS files only with `output_rds_path`. `create_html_report()` also writes
    `<output_file>_report_text_output.txt`; `create_pdfs_from_rds_objects()`
    writes the PDF files in the folder `pdf/` next to the RDS folder.
    `generate_rdata_list()` needs a named vector of paths (before: a list).
    The values of `get_metagene_data()` and `get_codon_usage_data()` are as
    the functions return them (before: other element names), and
    `get_metagene_data()` needs the `res_all` of a `_results_RiboseQC` file.
    `plot_frame_dist_boxplot()` plots the mean and its standard error (before:
    boxplots). `get_rl_and_cutoffs()` gives all read lengths with a cutoff
    (before: the selected ones), and `get_default_rl_selection()` the read
    lengths of `max_coverage`, which are not those of the P-sites with
    another method. `get_top50_all_genes()` gives 100 genes (before: 50). The
    codon usage is per codon position (before: per nucleotide position), and
    the genetic code comes from `RiboseQC_analysis()` (before: "a
    user-defined genetic code"). The links to `plot_codon_usage()`, which
    does not exist, now go to `plot_codon_usage_positional()`.
  - `?calc_cutoffs_from_profiles`: the value of the function is complete
    (before: an unescaped `%` cut the text), and it is `NULL` with fewer than
    5 transcripts. Four methods give the cutoffs (before: three), from 24 nt
    before the start codon to its first nt, so the cutoff is at most 24.
    `?choose_readlengths`: `summary_data` is `results_cutoffs` (before: the
    output of `calc_cutoffs_from_profiles()`), and `max_coverage` uses only
    the windows of the CDS.
- The new help page `?RiboseQC` gives the steps of an analysis and the
  functions for them (before, `?RiboseQC` found no page).
- The new help page `?RiboseQC_output` (also `?res_all`) describes each
  output file of `RiboseQC_analysis()`, each element of the list `res_all`
  (its rows, columns, units and the reads that it counts), the files of the
  report and the file of `prepare_for_ORFquant()`. The description of
  `res_all` moved there from `?RiboseQC_analysis`, and it now says what the
  code does where the old text did not: for example, `sequence_analysis`
  has the reads of a 5' end and a read length (the range has the width of
  the read length, and `seq` is the sequence of the genome), the P-sites of
  `_for_ORFquant` and of the bedgraph files are summed over the read
  lengths, so a range can be longer than 1 nt, `summary_P_sites` has all
  read lengths with an offset (not only the selected ones), and
  `P_sites_percodon_ratio` is `NaN` where a codon does not occur. The code
  does not change.

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

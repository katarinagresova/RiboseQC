# Input that cannot work stops RiboseQC_analysis() before it reads the annotation

test_that("RiboseQC_analysis() refuses wrong n_cores", {
  for (n in list(0, NA, "2", c(1, 2))) {
    expect_error(RiboseQC_analysis(annotation_file = "a", bam_files = "b", n_cores = n),
                 "n_cores must be a number")
  }
})

test_that("RiboseQC_analysis() refuses values of stranded other than TRUE, FALSE and 'inverse', before it reads a file", {
  # Before, "reverse" and "yes" stopped the run in the counting pass of the first BAM file with "'ignore.strand'
  # must be TRUE or FALSE", and NA with "missing value where TRUE/FALSE needed"
  for (s in list("reverse", "yes", NA, "Inverse")) {
    expect_error(RiboseQC_analysis(annotation_file = "a", bam_files = "b", stranded = s),
                 paste0("'stranded' must be TRUE, FALSE or 'inverse', not ", s), fixed = TRUE)
  }
  expect_error(RiboseQC_analysis(annotation_file = "a", bam_files = c("b1", "b2"), stranded = c(TRUE, "reverse")),
               "not reverse", fixed = TRUE)
  # the values that work go on to the next check: the files do not exist
  for (s in list(TRUE, FALSE, "inverse", "TRUE", "T", 1)) {
    expect_error(RiboseQC_analysis(annotation_file = "a", bam_files = "b", stranded = s), "don't exist")
  }
})

test_that("RiboseQC_analysis() refuses fewer annotation files than BAM files, but more than 1", {
  expect_error(RiboseQC_analysis(annotation_file = c("a1", "a2"), bam_files = c("b1", "b2", "b3")),
               "length of 'annotation_file' must be 1 or same as 'bam_files'")
})

test_that("RiboseQC_analysis() names all input files that do not exist, before the first BAM file, and writes nothing", {
  work <- withr::local_tempdir()
  msg <- tryCatch(RiboseQC_analysis(annotation_file = file.path(work, "no_Rannot"),
                                    bam_files = file.path(work, c("no1.bam", "no2.bam")),
                                    dest_names = file.path(work, c("s1", "s2")), create_report = FALSE),
                  error = conditionMessage)
  expect_match(msg, "don't exist")
  for (f in c("no_Rannot", "no1.bam", "no2.bam")) expect_match(msg, f, fixed = TRUE)
  # The first BAM file exists: before, the run stopped only after its analysis
  root <- system.file("ext_data", "simp_arab_root.bam", package = "RiboseQC")
  msg <- tryCatch(RiboseQC_analysis(annotation_file = root, bam_files = c(root, file.path(work, "no2.bam")),
                                    dest_names = file.path(work, c("s1", "s2")), create_report = FALSE),
                  error = conditionMessage)
  expect_match(msg, "don't exist:\n[^\n]*no2\\.bam\n$")
  expect_length(list.files(work, all.files = TRUE, no.. = TRUE), 0)
})

test_that("RiboseQC_analysis() stops before the first BAM file when the FASTA file of the annotation is gone", {
  work <- withr::local_tempdir()
  ext <- system.file("ext_data", package = "RiboseQC")
  fasta <- file.path(work, "fasta", "simp_arab.fa.gz")
  dir.create(dirname(fasta))
  file.copy(file.path(ext, paste0("simp_arab.fa.gz", c("", ".fai", ".gzi"))), dirname(fasta))
  quiet(annotation <- prepare_annotation_files(annotation_directory = file.path(work, "annotation"),
                                               gtf_file = file.path(ext, "simp_arab.gtf.gz"), genome_seq = fasta),
        expected = c("some CDS phases are missing", "genome version information is not available"))
  unlink(dirname(fasta), recursive = TRUE)
  # Before, the run stopped in the P-sites profiles, near the end of the analysis of the first BAM file, with
  # "file(s) do not exist"
  msg <- tryCatch(suppressMessages(utils::capture.output(
    RiboseQC_analysis(annotation_file = annotation, bam_files = file.path(ext, "simp_arab_root.bam"),
                      dest_names = file.path(work, "out", "root"), create_report = FALSE))),
    error = conditionMessage)
  expect_match(msg, "The FASTA file of the annotation .*simp_arab.gtf.gz_Rannot doesn't exist: .*fasta/simp_arab\\.fa\\.gz")
  expect_false(dir.exists(file.path(work, "out")))
})

# The annotation of the example data, made once for the tests below
ext <- system.file("ext_data", package = "RiboseQC")
bams <- file.path(ext, c("simp_arab_root.bam", "simp_arab_shoots.bam"))
work <- tempfile("riboseqc_analysis_input_")
dir.create(work)
example_annotation <- local({
  annotation <- NULL
  function() {
    if (is.null(annotation)) {
      quiet(annotation <<- prepare_annotation_files(annotation_directory = file.path(work, "annotation"),
                                                    gtf_file = file.path(ext, "simp_arab.gtf.gz"),
                                                    genome_seq = file.path(ext, "simp_arab.fa.gz")),
            expected = c("some CDS phases are missing", "genome version information is not available"))
    }
    annotation
  }
})
# A copy of the annotation with the seqinfo change(seqinfo); its path
changed_annotation <- function(change) {
  ann <- get(load(example_annotation()))
  ann$seqinfo <- change(ann$seqinfo)
  f <- tempfile("changed_Rannot_", work)
  save(ann, file = f)
  f
}
# The error of RiboseQC_analysis(), without its output
analysis_error <- function(..., create_report = FALSE) {
  tryCatch({
    suppressMessages(utils::capture.output(RiboseQC_analysis(..., create_report = create_report)))
    NA_character_
  }, error = conditionMessage)
}

test_that("RiboseQC_analysis() refuses a BAM file with other chromosomes than the annotation, before the first BAM file", {
  other <- changed_annotation(function(si) {
    Seqinfo::seqlevels(si) <- paste0("other", seq_along(Seqinfo::seqlevels(si)))
    si
  })
  out <- file.path(work, "out_names")
  # Before, a BAM file with other chromosome names than its annotation stopped the run in the analysis of the BAM
  # file, with "1 elements in value to replace 0 elements". This annotation, with other names in its seqinfo only,
  # gave no error
  msg <- analysis_error(annotation_file = other, bam_files = bams[1], dest_names = file.path(out, "root"))
  expect_match(msg, "No chromosome of the BAM file is in the annotation")
  expect_match(msg, "simp_arab_root.bam has e.g. ChrC, ChrM, Chr4", fixed = TRUE)
  expect_match(msg, "not to its transcripts")
  # The second BAM file has its own annotation, which has other chromosomes
  msg <- analysis_error(annotation_file = c(example_annotation(), other), bam_files = bams,
                        dest_names = file.path(out, c("root", "shoots")))
  expect_match(msg, "the BAM file .*simp_arab_shoots.bam has")
  expect_false(dir.exists(out))
})

test_that("RiboseQC_analysis() refuses a BAM file with other chromosome lengths than the annotation", {
  longer <- changed_annotation(function(si) {
    Seqinfo::seqlengths(si) <- Seqinfo::seqlengths(si) + 1L
    si
  })
  out <- file.path(work, "out_lengths")
  # Before, a BAM file of another genome assembly stopped the run in the analysis of the BAM file, with
  # "incompatible seqlengths". This annotation, with other lengths in its seqinfo only, gave no error
  msg <- analysis_error(annotation_file = longer, bam_files = bams[1], dest_names = file.path(out, "root"))
  expect_match(msg, "is not on the genome of the annotation")
  expect_match(msg, "7 chromosomes have other lengths, e.g. ChrC has 154478 nt in the BAM file and 154479 nt in the annotation",
               fixed = TRUE)
  expect_false(dir.exists(out))
})

test_that("RiboseQC_analysis() does not use genome_seq, and warns when it is not the FASTA file of the annotation", {
  # with other chromosomes in the annotation, the run stops after the annotation is loaded
  other <- changed_annotation(function(si) {
    Seqinfo::seqlevels(si) <- paste0("other", seq_along(Seqinfo::seqlevels(si)))
    si
  })
  run <- function(genome_seq) {
    warns <- character()
    msg <- withCallingHandlers(
      analysis_error(annotation_file = other, bam_files = bams[1], dest_names = file.path(work, "out_genome", "root"),
                     genome_seq = genome_seq),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    list(msg = msg, warns = grep("genome_seq", warns, value = TRUE))
  }
  fasta <- file.path(ext, "simp_arab.fa.gz")
  # Before, an FaFile stopped the run with "invalid 'file' argument"
  for (genome_seq in list(NULL, fasta, Rsamtools::FaFile(fasta))) {
    res <- run(genome_seq)
    expect_match(res$msg, "No chromosome of the BAM file is in the annotation")
    expect_length(res$warns, 0)
  }
  other_fasta <- file.path(work, "other.fa")
  for (genome_seq in list(other_fasta, Rsamtools::FaFile(other_fasta))) {
    res <- run(genome_seq)
    expect_match(res$msg, "No chromosome of the BAM file is in the annotation")
    expect_length(res$warns, 1)
    expect_match(res$warns, "genome_seq \\(.*other\\.fa\\) is not used: .*the FASTA file .*simp_arab\\.fa\\.gz")
  }
  # a relative path to the FASTA file of the annotation
  withr::with_dir(ext, expect_silent(check_genome_seq("simp_arab.fa.gz", Rsamtools::FaFile(fasta), "annotation")))
})

test_that("RiboseQC_analysis() refuses output directories that cannot be written to, and makes those that do not exist", {
  skip_on_os("windows")
  read_only <- file.path(work, "read_only")
  dir.create(read_only)
  Sys.chmod(read_only, "555")
  withr::defer(Sys.chmod(read_only, "755"))
  skip_if(file.access(read_only, 2) == 0, "the user can write to a read-only directory")
  # Before, the run stopped at the first file that it wrote, with "cannot open the connection"
  msg <- analysis_error(annotation_file = example_annotation(), bam_files = bams[1],
                        dest_names = file.path(read_only, "root"))
  expect_match(msg, "Cannot write the output files, .*read_only is not a directory that exists and can be written to")
  # The directory of the report. Before, the run stopped at the report, after the analysis. The output
  # directory that does not exist is made, as before
  out <- file.path(work, "out_report", "sub")
  msg <- analysis_error(annotation_file = example_annotation(), bam_files = bams[1], dest_names = file.path(out, "root"),
                        create_report = TRUE, report_file = file.path(read_only, "report.html"))
  expect_match(msg, "Cannot write the output files, .*read_only is not a directory .*report\\.html")
  expect_true(dir.exists(out))
  expect_length(list.files(out, all.files = TRUE, no.. = TRUE), 0)
})

test_that("RiboseQC_analysis() says when a BAM file has no mapped reads", {
  empty <- file.path(work, "empty.bam")
  Rsamtools::filterBam(bams[1], empty, filter = S4Vectors::FilterRules(list(none = function(x) rep(FALSE, nrow(x)))))
  # Before, the run stopped with "subscript out of bounds"
  msg <- analysis_error(annotation_file = example_annotation(), bam_files = empty,
                        dest_names = file.path(work, "out_empty", "empty"))
  expect_match(msg, "The BAM file .*empty\\.bam has no mapped reads")
})

test_that("RiboseQC_analysis() warns for the reads on chromosomes that the annotation does not have", {
  si <- Seqinfo::seqinfo(Rsamtools::BamFile(bams[1]))
  expect_silent(warn_bam_chromosomes(bams[1], si, "annot"))
  # The root BAM file has 1895 alignments on ChrC and 6787 on ChrM
  expect_warning(warn_bam_chromosomes(bams[1], si[setdiff(Seqinfo::seqlevels(si), c("ChrC", "ChrM"))], "annot"),
                 paste("simp_arab_root.bam has 8682 alignments on 2 chromosomes that the annotation annot does not have,",
                       "e.g. ChrC (1895), ChrM (6787). RiboseQC_analysis() uses these reads"), fixed = TRUE)
  # The other chromosomes of this BAM file have no reads
  chr1 <- file.path(work, "chr1.bam")
  Rsamtools::filterBam(bams[1], chr1, param = Rsamtools::ScanBamParam(which = GenomicRanges::GRanges("Chr1", IRanges::IRanges(1, 30427671))))
  expect_silent(warn_bam_chromosomes(chr1, si["Chr1"], "annot"))
  # A BAM file without an index: the chromosomes of its header. RiboseQC_analysis() warns before the first BAM
  # file; this one has no reads, so the run stops at the start of its analysis
  empty <- file.path(work, "empty_noindex.bam")
  Rsamtools::filterBam(bams[1], empty, indexDestination = FALSE,
                       filter = S4Vectors::FilterRules(list(none = function(x) rep(FALSE, nrow(x)))))
  no_chrc <- changed_annotation(function(si) si[setdiff(Seqinfo::seqlevels(si), "ChrC")])
  warns <- character()
  msg <- withCallingHandlers(
    analysis_error(annotation_file = no_chrc, bam_files = empty, dest_names = file.path(work, "out_chrc", "empty")),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_match(msg, "empty_noindex\\.bam has no mapped reads")
  expect_length(warns, 1)
  expect_match(warns, "empty_noindex\\.bam has 1 chromosome that the annotation .*changed_Rannot_.* does not have, e\\.g\\. ChrC\\. The BAM file has no index")
})

# prepare_for_ORFquant() on the example data (7 genes of human chr22, see
# data/README.md). From ORFquant's tests (tests/testthat/test-pipeline.R and
# test-input-checks.R), without run_ORFquant().

skip_on_os("windows")   # with n_cores > 1, prepare_for_ORFquant() forks

work <- tempfile("riboseqc_for_orfquant_")
dir.create(work)
runs <- new.env()   # the annotation file, and the prefix of the files of each run
quiet <- function(expr) invisible(utils::capture.output(suppressMessages(suppressWarnings(expr))))
# prepare_for_ORFquant() leaves its annotation in two global variables
withr::defer(rm(list = intersect(c("GTF_annotation", "genome_seq"), ls(globalenv())),
                envir = globalenv()), teardown_env())

test_that("prepare_for_ORFquant() refuses wrong n_cores", {
  expect_error(prepare_for_ORFquant(annotation_file = "a", bam_file = "b", path_to_rl_cutoff_file = "c",
                                    n_cores = 0), "n_cores must be a number")
})

test_that("prepare_for_ORFquant() names the files that do not exist", {
  msg <- tryCatch(prepare_for_ORFquant(annotation_file = file.path(work, "no_Rannot"),
                                       bam_file = file.path(work, "no.bam"),
                                       path_to_rl_cutoff_file = file.path(work, "no.tsv")),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "don't exist")
  for (f in c("no_Rannot", "no.bam", "no.tsv")) expect_match(msg, f, fixed = TRUE)
})

# Calculates the P-sites of the example; the prefix
prepare_example <- function(name, n_cores, chunk_size) {
  prefix <- file.path(work, name)
  quiet(prepare_for_ORFquant(annotation_file = runs$annotation,
                             bam_file = example_file("chr22_example.bam"),
                             path_to_rl_cutoff_file = example_file("chr22_example_cutoffs.tsv"),
                             dest_name = prefix, n_cores = n_cores, chunk_size = chunk_size))
  prefix
}

test_that("the annotation of the example is prepared", {
  quiet(prepare_annotation_files(annotation_directory = work,
                                 gtf_file = example_file("chr22_example.gtf.gz"),
                                 genome_seq = example_file("chr22_example.fa.gz"),
                                 forge_BSgenome = FALSE))
  runs$annotation <- file.path(work, "chr22_example.gtf.gz_Rannot")
  expect_true(file.exists(runs$annotation))
})

test_that("prepare_for_ORFquant() gives the P-sites and the junction reads of the example", {
  skip_if(is.null(runs$annotation))
  rm(list = intersect(c("GTF_annotation", "genome_seq"), ls(globalenv())), envir = globalenv())
  runs$one_core <- prepare_example("one_core", n_cores = 1, chunk_size = 5e6)
  res <- get(load(paste0(runs$one_core, "_for_ORFquant")))
  expect_named(res, c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm", "junctions"))
  # The numbers of positions and of P-sites, and of junction reads, that prepare_for_ORFquant()
  # gave in ORFquant, before it moved to RiboseQC
  sizes <- vapply(res[c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm")],
                  function(g) c(length(g), sum(g$score)), numeric(2))
  expect_equal(unname(sizes), matrix(c(370, 506, 302, 406, 26, 27), nrow = 2))
  expect_equal(c(length(res$junctions), sum(res$junctions$reads), sum(res$junctions$unique_reads)),
               c(186, 101, 77))
  # it leaves the annotation in the global variables, as before
  expect_true(all(c("GTF_annotation", "genome_seq") %in% ls(globalenv())))
})

test_that("the P-sites do not depend on the number of cores or of chunks", {
  skip_if(is.null(runs$one_core))
  # 525 alignments in the BAM file: the P-sites of 6 chunks are calculated
  # in forked processes
  runs$two_cores <- prepare_example("two_cores", n_cores = 2, chunk_size = 100)
  expect_identical(get(load(paste0(runs$one_core, "_for_ORFquant"))),
                   get(load(paste0(runs$two_cores, "_for_ORFquant"))))
})

# Input that cannot work stops the run before the BAM file is read

# The error message of a call, and the warnings that came before it
refusal <- function(expr) {
  warned <- character()
  msg <- tryCatch(withCallingHandlers(suppressMessages(expr),
                                      warning = function(w) {
                                        warned <<- c(warned, conditionMessage(w))
                                        invokeRestart("muffleWarning")
                                      }),
                  error = function(e) conditionMessage(e))
  list(message = if (is.character(msg)) msg else NA_character_, warnings = warned)
}
# The example annotation with changed seqinfo
changed_annotation <- function(change) {
  ann <- get(load(runs$annotation))
  ann$seqinfo <- change(ann$seqinfo)
  f <- tempfile("changed_Rannot_", work)
  save(ann, file = f)
  f
}
prepare_with <- function(cutoffs, annotation = runs$annotation, dest = file.path(work, "refused")) {
  refusal(prepare_for_ORFquant(annotation_file = annotation,
                               bam_file = example_file("chr22_example.bam"),
                               path_to_rl_cutoff_file = cutoffs, dest_name = dest))
}
write_cutoffs <- function(lines) {
  f <- tempfile("cutoffs_", work, fileext = ".tsv")
  writeLines(lines, f)
  f
}

test_that("prepare_for_ORFquant() refuses an output directory that does not exist", {
  skip_if(is.null(runs$annotation))
  msg <- prepare_with(example_file("chr22_example_cutoffs.tsv"), dest = file.path(work, "no_dir", "sample"))$message
  expect_match(msg, "Cannot write the output files")
})

test_that("prepare_for_ORFquant() refuses a wrong cutoff table", {
  skip_if(is.null(runs$annotation))
  head <- "read_length\tcutoff\tcompartment"
  expect_match(prepare_with(write_cutoffs(head))$message, "has no rows")
  expect_match(prepare_with(write_cutoffs(c(head, "twenty\t12\tnucl")))$message, "must be numbers")
  out <- prepare_with(write_cutoffs(c(head, "28\t12\tmito")))
  expect_match(out$message, "No compartment of the rl_cutoff file is 'nucl'")
  expect_match(out$message, "mito")
})

test_that("prepare_for_ORFquant() refuses a BAM file with other chromosomes than the annotation", {
  skip_if(is.null(runs$annotation))
  other <- changed_annotation(function(si) {
    Seqinfo::seqlevels(si) <- paste0("other", seq_along(Seqinfo::seqlevels(si)))
    si
  })
  out <- prepare_with(write_cutoffs(c("read_length\tcutoff\tcompartment", "28\t12\tnucl", "29\t12\tmito")),
                      annotation = other)
  expect_match(out$message, "No chromosome of the BAM file is in the annotation")
  expect_match(out$message, "chr22")
  # the compartment that is not used is a warning, not an error
  expect_match(out$warnings, "mito")
})

test_that("prepare_for_ORFquant() refuses a BAM file with other chromosome lengths than the annotation", {
  skip_if(is.null(runs$annotation))
  # Before, it stopped in the first chunk of reads with "sequence chr22 has incompatible seqlengths"
  longer <- changed_annotation(function(si) {
    Seqinfo::seqlengths(si) <- Seqinfo::seqlengths(si) + 1L
    si
  })
  msg <- prepare_with(example_file("chr22_example_cutoffs.tsv"), annotation = longer)$message
  expect_match(msg, "is not on the genome of the annotation")
  expect_match(msg, "chr22 has [0-9]+ nt in the BAM file and [0-9]+ nt in the annotation")
})

test_that("prepare_for_ORFquant() refuses a cutoff table that gives no P-sites", {
  skip_if(is.null(runs$annotation))
  msg <- prepare_with(write_cutoffs(c("read_length\tcutoff\tcompartment", "99\t12\tnucl")))$message
  expect_match(msg, "No P-sites were found")
  expect_false(file.exists(file.path(work, "refused_for_ORFquant")))
})

test_that("prepare_for_ORFquant() sorts the P-sites of bigWig files by the chromosomes of the annotation", {
  skip_if(is.null(runs$one_core))
  # The annotation has chr22 and then chr1; the bigWig files have chr1 first
  annotation <- changed_annotation(function(si) {
    Seqinfo::seqlevels(si) <- c(Seqinfo::seqlevels(si), "chr1")
    Seqinfo::seqlengths(si)[["chr1"]] <- Seqinfo::seqlengths(si)[["chr22"]]
    si
  })
  ps <- get(load(paste0(runs$one_core, "_for_ORFquant")))$P_sites_all
  in_bw <- GenomicRanges::GRanges(rep(c("chr1", "chr22"), each = length(ps)), rep(IRanges::ranges(ps), 2),
                                  rep(BiocGenerics::strand(ps), 2), score = rep(ps$score, 2),
                                  seqinfo = get(load(annotation))$seqinfo[c("chr1", "chr22")])
  bw <- file.path(work, c("plus.bw", "minus.bw"))
  rtracklayer::export(in_bw[BiocGenerics::strand(in_bw) == "+"], bw[1], format = "bigWig")
  rtracklayer::export(in_bw[BiocGenerics::strand(in_bw) == "-"], bw[2], format = "bigWig")
  dest <- file.path(work, "bigwig")
  quiet(prepare_for_ORFquant(annotation_file = annotation, bam_file = example_file("chr22_example.bam"),
                             path_to_P_sites_plus_bw = bw[1], path_to_P_sites_minus_bw = bw[2], dest_name = dest))
  got <- get(load(paste0(dest, "_for_ORFquant")))$P_sites_all
  expect_false(BiocGenerics::is.unsorted(got))
  expected <- in_bw
  Seqinfo::seqlevels(expected) <- c("chr22", "chr1")
  expected <- BiocGenerics::sort(expected)
  expected$score <- as.numeric(expected$score)   # the scores of bigWig files are double
  expect_identical(as.data.frame(got), as.data.frame(expected))
})

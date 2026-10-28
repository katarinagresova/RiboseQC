# The command-line script inst/scripts/run_riboseqc.R: what it checks before it starts
# the steps, and a run of the steps on the example data

cli_script <- system.file("scripts", "run_riboseqc.R", package = "RiboseQC")
rscript <- file.path(R.home("bin"), "Rscript")
ext_data <- system.file("ext_data", package = "RiboseQC")

# The exit status and the output (standard output and error) of the script
run_cli <- function(...) {
  out <- suppressWarnings(system2(rscript, shQuote(c(cli_script, ...)), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = paste(out, collapse = "\n"))
}

# The options of a run that is correct, as a vector of arguments. An argument of this
# function replaces an option, and NULL removes it.
cli_args <- function(...) {
  opts <- list(gtf = file.path(ext_data, "simp_arab.gtf.gz"),
               fasta = file.path(ext_data, "simp_arab.fa.gz"),
               bam = file.path(ext_data, "simp_arab_root.bam"),
               outdir = tempfile("riboseqc_cli_"))
  opts <- utils::modifyList(opts, list(...))
  as.vector(rbind(paste0("--", names(opts)), unlist(opts)))
}

test_that("--help prints the usage and succeeds", {
  r <- run_cli("--help")
  expect_equal(r$status, 0L)
  expect_match(r$output, "^Usage: Rscript run_riboseqc.R")
  expect_match(r$output, "--annotation FILE", fixed = TRUE)
})

test_that("without arguments it prints the usage and fails", {
  r <- run_cli()
  expect_equal(r$status, 1L)
  expect_match(r$output, "^Usage: Rscript run_riboseqc.R")
})

test_that("it stops with a message on wrong options, before it runs a step", {
  expect_cli_error <- function(args, message) {
    r <- run_cli(args)
    expect_equal(r$status, 1L, info = r$output)
    expect_match(r$output, message, fixed = TRUE)
    expect_false(grepl("failed", r$output, fixed = TRUE))   # no step ran
  }
  gtf <- file.path(ext_data, "simp_arab.gtf.gz")
  bam <- file.path(ext_data, "simp_arab_root.bam")
  expect_cli_error(c(cli_args(), "--nonsense", "1"), "unknown option --nonsense")
  expect_cli_error(cli_args(fasta = NULL), "missing --fasta")
  expect_cli_error(cli_args(bam = NULL), "missing --bam")
  expect_cli_error(cli_args(annotation = gtf), "give either --gtf and --fasta, or --annotation")
  expect_cli_error(cli_args(gtf = NULL, fasta = NULL), "give either --gtf and --fasta, or --annotation")
  expect_cli_error(cli_args(bam = paste0(bam, ",no_such_file.bam")), "--bam: file not found: no_such_file.bam")
  expect_cli_error(cli_args(bam = paste0(bam, ",", bam)), "two BAM files have the sample name simp_arab_root")
  expect_cli_error(cli_args("sample-names" = "a,b"),
                   "--sample-names must have one name for each of the 1 BAM files, not 2")
  expect_cli_error(cli_args(stranded = "yes"), "--stranded must be TRUE, FALSE or inverse, not yes")
  expect_cli_error(cli_args("readlength-method" = "best"),
                   "--readlength-method must be max_coverage, max_inframe or all, not best")
  expect_cli_error(cli_args(cores = "0"), "--cores must be a whole number above 0, not 0")
  expect_cli_error(cli_args(cores = "two"), "--cores must be a whole number above 0, not two")
  expect_cli_error(c(cli_args(), "--no-report", "--extended-report"), "give either --no-report or --extended-report")
  expect_cli_error(c(cli_args(), "--no-report=yes"), "--no-report takes no value")
})

test_that("the command-line script runs with the annotation, and gives the results of RiboseQC_analysis()", {
  # Its steps run in new R processes, which need RiboseQC to be installed
  loads <- system2(rscript, c("-e", shQuote("quit(status = !requireNamespace('RiboseQC', quietly = TRUE))")),
                   stdout = FALSE, stderr = FALSE) == 0
  skip_if_not(loads, "RiboseQC is not installed in the library of Rscript")
  example <- example_analysis()
  outdir <- tempfile("riboseqc_cli_")
  # the shoots BAM file: the root sample has no P-sites, so no _for_ORFquant
  r <- run_cli("--annotation", example$annotation, "--bam", file.path(ext_data, "simp_arab_shoots.bam"),
               "--sample-names", "shoots", "--outdir", outdir, "--cores", "2", "--no-report")
  expect_equal(r$status, 0L, info = r$output)
  # the list of the output files at the end
  listed <- read.delim(text = sub(".*RiboseQC output files:\n", "", r$output), header = FALSE)
  expect_equal(listed$V2[listed$V1 %in% c("annotation", "offsets", "for_ORFquant")],
               c(normalizePath(example$annotation),
                 file.path(normalizePath(outdir), c("shoots_P_sites_offsets.tsv", "shoots_for_ORFquant"))))
  expect_false(any(c("report", "report_plots") %in% listed$V1))
  # all tab-separated files, with those of write_TSV_files = TRUE (the default)
  expect_setequal(listed$V1[grepl("\\.tsv$", listed$V2)],
                  c("offsets", "summary", "read_lengths", "regions", "genes", "P_sites_calcs_tsv", "top50_positions"))
  expect_setequal(listed$V2[grepl("\\.tsv$", listed$V2)],
                  file.path(normalizePath(outdir), list.files(outdir, pattern = "\\.tsv$")))
  # the same files as the run of RiboseQC_analysis() in helper-example.R
  for (f in c("_P_sites_offsets.tsv", "_P_sites_calcs", "_P_sites_plus.bedgraph", "_coverage_plus.bedgraph")) {
    expect_identical(readLines(file.path(outdir, paste0("shoots", f))), readLines(paste0(example$prefix[["shoots"]], f)),
                     info = f)
  }
  # equal, not identical: the sequences of the top positions (a DNAStringSet) are in a memory pool that
  # is different in each run
  expect_equal(get(load(file.path(outdir, "shoots_results_RiboseQC_all"))), example_results("shoots"), tolerance = 0)
})

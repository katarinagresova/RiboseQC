# A file of the test data in tests/testthat/data (see its README.md)
example_file <- function(f) test_path("data", f)

# Runs `expr` without its output and messages. A warning fails the test, unless it
# contains one of the texts `expected`
quiet <- function(expr, expected = character()) {
  warned <- character()
  utils::capture.output(suppressMessages(withCallingHandlers(expr, warning = function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  })))
  for (e in expected) warned <- warned[!grepl(e, warned, fixed = TRUE)]
  expect(length(warned) == 0, paste(c("unexpected warnings:", warned), collapse = "\n"))
}

# The example data of the package (inst/ext_data): its annotation, and
# RiboseQC_analysis() with the default parameters, without the report, on its two BAM
# files. They run once, at the first call, for all tests (about 3 min). The path of the
# annotation, and the prefix of the output files of each sample
example_analysis <- local({
  res <- NULL
  function() {
    if (is.null(res)) {
      ext_data <- system.file("ext_data", package = "RiboseQC")
      work <- tempfile("riboseqc_example_")
      dir.create(work)
      # txdbmaker warns that some CDS phases are missing, and that the TxDb has no
      # genome version
      quiet(annotation <- prepare_annotation_files(annotation_directory = work,
                                                   gtf_file = file.path(ext_data, "simp_arab.gtf.gz"),
                                                   genome_seq = file.path(ext_data, "simp_arab.fa.gz")),
            expected = c("some CDS phases are missing", "genome version information is not available"))
      prefix <- c(root = file.path(work, "root"), shoots = file.path(work, "shoots"))
      quiet(RiboseQC_analysis(annotation_file = annotation,
                              bam_files = file.path(ext_data, c("simp_arab_root.bam", "simp_arab_shoots.bam")),
                              dest_names = prefix, create_report = FALSE))
      res <<- list(annotation = annotation, prefix = prefix)
    }
    res
  }
})

# The results of RiboseQC_analysis() on `sample` of the example data: res_all of its
# _results_RiboseQC_all file
example_results <- function(sample) {
  get(load(paste0(example_analysis()$prefix[[sample]], "_results_RiboseQC_all")))
}

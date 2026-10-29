# RiboseQC_analysis() gives the same results in 2 chunks as in one chunk. A chunk has more
# than 100,000 alignments (chunk_size), and the BAM files of the example data have fewer:
# so the test uses 2 copies of the shoots BAM file (106,284 alignments). Chr1 is the last
# chromosome of the BAM file: the second chunk has 6,283 alignments on Chr1, and the two
# copies of the last read of the first chunk are in different chunks. With
# rescue_all_rls = TRUE, the read lengths of nucl get the offset 12, so both chunks have
# P-sites. About 3 min.

test_that("RiboseQC_analysis() gives the same files in 2 chunks as in one chunk", {
  work <- tempfile("riboseqc_chunks_")
  dir.create(work)
  bam <- file.path(work, "shoots_twice.bam")
  Rsamtools::mergeBam(rep(system.file("ext_data", "simp_arab_shoots.bam", package = "RiboseQC"), 2), bam,
                      indexDestination = TRUE)
  expect_equal(Rsamtools::countBam(bam)$records, 106284)
  prefix <- c(one = file.path(work, "one"), two = file.path(work, "two"))
  chunk_size <- c(one = 5e6, two = 100001)
  for (n in names(prefix)) {
    quiet(RiboseQC_analysis(annotation_file = example_analysis()$annotation, bam_files = bam,
                            dest_names = prefix[[n]], chunk_size = chunk_size[[n]], rescue_all_rls = TRUE,
                            create_report = FALSE))
  }
  files <- sub("^one", "", list.files(work, "^one_"))
  expect_setequal(sub("^two", "", list.files(work, "^two_")), files)
  for (f in files) {
    if (f %in% c("_results_RiboseQC", "_results_RiboseQC_all", "_for_ORFquant", "_junctions")) {
      # equal, not identical: the sequences of the top positions (a DNAStringSet) are in a memory
      # pool that is different in each run
      expect_equal(get(load(paste0(prefix[["two"]], f))), get(load(paste0(prefix[["one"]], f))),
                   tolerance = 0, info = f)
    } else {
      expect_identical(readLines(paste0(prefix[["two"]], f)), readLines(paste0(prefix[["one"]], f)), info = f)
    }
  }
})

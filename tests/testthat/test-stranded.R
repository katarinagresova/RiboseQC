# RiboseQC_analysis(stranded = "inverse") on the root BAM file of the example data
# with the strand of each read inverted must give the results of the root BAM file
# with stranded = TRUE (example_analysis() in helper-example.R). The run also has
# write_TSV_files = FALSE, which does not change the results

test_that("RiboseQC_analysis(stranded = \"inverse\") gives the results of stranded = TRUE on the reads with inverted strands", {
  work <- withr::local_tempdir()
  inverted <- invert_bam(system.file("ext_data", "simp_arab_root.bam", package = "RiboseQC"),
                         file.path(work, "root_inverted.bam"))
  quiet(RiboseQC_analysis(annotation_file = example_analysis()$annotation, bam_files = inverted,
                          dest_names = file.path(work, "root"), create_report = FALSE, stranded = "inverse",
                          write_TSV_files = FALSE))
  # Before, the P-sites, the coverage and the junctions were on the strands of the BAM
  # file, and the reads of each region and gene were counted on both strands
  expect_equal(get(load(file.path(work, "root_results_RiboseQC_all"))), example_results("root"))
  # of the tab-separated files, only P_sites_offsets.tsv and P_sites_calcs, which do not depend on write_TSV_files
  expect_setequal(list.files(work, pattern = "tsv$|calcs"), c("root_P_sites_offsets.tsv", "root_P_sites_calcs"))
})

# With stranded = FALSE, the reads of each region by read length ignore the strand, as the
# reads of each region (positions). So the root BAM file with the strand of each read inverted
# gives the counts of the root BAM file with stranded = TRUE: the bundled reads are not on
# regions of the other strand
test_that("RiboseQC_analysis(stranded = FALSE) counts the reads of each region by read length on both strands", {
  work <- withr::local_tempdir()
  inverted <- invert_bam(system.file("ext_data", "simp_arab_root.bam", package = "RiboseQC"),
                         file.path(work, "root_inverted.bam"))
  # The 5' ends of the inverted reads are on the other strand of the transcripts, so they have no frame
  quiet(RiboseQC_analysis(annotation_file = example_analysis()$annotation, bam_files = inverted,
                          dest_names = file.path(work, "root"), create_report = FALSE, stranded = FALSE,
                          write_TSV_files = FALSE),
        expected = "Not enough signal or low frame preference, skipped P-sites")
  rs <- get(load(file.path(work, "root_results_RiboseQC_all")))$read_stats
  expected <- example_results("root")$read_stats
  expect_equal(rs$positions, expected$positions)
  # Before, reads_summary counted the reads only on the regions of their strand
  expect_equal(rs$reads_summary, expected$reads_summary)
  expect_equal(rs$reads_summary_unq, expected$reads_summary_unq)
})

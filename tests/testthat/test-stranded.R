# RiboseQC_analysis(stranded = "inverse") on the shoots BAM file of the example data
# with the strand of each read inverted must give the results of the shoots BAM file
# with stranded = TRUE (example_analysis() in helper-example.R). The run also has
# write_TSV_files = FALSE, which does not change the results. The shoots sample has
# P-sites (on ChrC); the root sample has none

test_that("RiboseQC_analysis(stranded = \"inverse\") gives the results of stranded = TRUE on the reads with inverted strands", {
  work <- withr::local_tempdir()
  inverted <- invert_bam(system.file("ext_data", "simp_arab_shoots.bam", package = "RiboseQC"),
                         file.path(work, "shoots_inverted.bam"))
  quiet(RiboseQC_analysis(annotation_file = example_analysis()$annotation, bam_files = inverted,
                          dest_names = file.path(work, "shoots"), create_report = FALSE, stranded = "inverse",
                          write_TSV_files = FALSE))
  # Before, the P-sites, the coverage and the junctions were on the strands of the BAM
  # file, and the reads of each region and gene were counted on both strands
  expect_equal(get(load(file.path(work, "shoots_results_RiboseQC_all"))), example_results("shoots"))
  # of the tab-separated files, only P_sites_offsets.tsv and P_sites_calcs, which do not depend on write_TSV_files
  expect_setequal(list.files(work, pattern = "tsv$|calcs"), c("shoots_P_sites_offsets.tsv", "shoots_P_sites_calcs"))
})

# RiboseQC_analysis(stranded = "inverse") on the root BAM file of the example data
# with the strand of each read inverted must give the results of the root BAM file
# with stranded = TRUE (example_analysis() in helper-example.R)

test_that("RiboseQC_analysis(stranded = \"inverse\") gives the results of stranded = TRUE on the reads with inverted strands", {
  work <- withr::local_tempdir()
  inverted <- invert_bam(system.file("ext_data", "simp_arab_root.bam", package = "RiboseQC"),
                         file.path(work, "root_inverted.bam"))
  quiet(RiboseQC_analysis(annotation_file = example_analysis()$annotation, bam_files = inverted,
                          dest_names = file.path(work, "root"), create_report = FALSE, stranded = "inverse"))
  # Before, the P-sites, the coverage and the junctions were on the strands of the BAM
  # file, and the reads of each region and gene were counted on both strands
  expect_equal(get(load(file.path(work, "root_results_RiboseQC_all"))), example_results("root"))
})

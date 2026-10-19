# Input that cannot work stops RiboseQC_analysis() before it reads the annotation

test_that("RiboseQC_analysis() refuses wrong n_cores", {
  for (n in list(0, NA, "2", c(1, 2))) {
    expect_error(RiboseQC_analysis(annotation_file = "a", bam_files = "b", n_cores = n),
                 "n_cores must be a number")
  }
})

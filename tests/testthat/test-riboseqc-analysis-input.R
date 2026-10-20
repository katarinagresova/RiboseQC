# Input that cannot work stops RiboseQC_analysis() before it reads the annotation

test_that("RiboseQC_analysis() refuses wrong n_cores", {
  for (n in list(0, NA, "2", c(1, 2))) {
    expect_error(RiboseQC_analysis(annotation_file = "a", bam_files = "b", n_cores = n),
                 "n_cores must be a number")
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

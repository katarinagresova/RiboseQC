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

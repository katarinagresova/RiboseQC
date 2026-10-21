# RiboseQC_analysis(stranded = "inverse") on the root BAM file of the example data
# with the strand of each read inverted must give the results of the root BAM file
# with stranded = TRUE (example_analysis() in helper-example.R)

# Writes the alignments of `bam` with the strand of each read inverted to the BAM
# file `dest`, and gives its path. All reads of the example data are mapped.
# RiboseQC reads the flag, the position, the CIGAR, MAPQ and the MD tag; the copy
# has no read names, sequences or qualities.
invert_bam <- function(bam, dest) {
  aln <- Rsamtools::scanBam(bam, param = Rsamtools::ScanBamParam(what = c("flag", "rname", "pos", "mapq", "cigar"),
                                                                 tag = "MD"))[[1]]
  chroms <- Rsamtools::scanBamHeader(bam)[[1]]$targets
  # scanBam() gives NA for MAPQ 255 (not available)
  mapq <- ifelse(is.na(aln$mapq), 255L, aln$mapq)
  md <- ifelse(is.na(aln$tag$MD), "", paste0("\tMD:Z:", aln$tag$MD))
  sam <- tempfile(fileext = ".sam")
  writeLines(c(paste0("@SQ\tSN:", names(chroms), "\tLN:", chroms),
               paste0("*\t", bitwXor(aln$flag, 16L), "\t", aln$rname, "\t", aln$pos, "\t", mapq, "\t", aln$cigar,
                      "\t*\t0\t0\t*\t*", md)), sam)
  Rsamtools::asBam(sam, sub("\\.bam$", "", dest))
}

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

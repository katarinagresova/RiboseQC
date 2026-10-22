# FaFile_Circ: a FASTA file in which ranges can go across the end of a circular
# chromosome. A genome of chrM (circular, 10 nt) and chr2 (not circular)

fasta_circ <- function() {
  f <- withr::local_tempfile(fileext = ".fa", .local_envir = parent.frame())
  Biostrings::writeXStringSet(setNames(DNAStringSet(c("AAAAAAAAGG", "AAAAAAAAGG")), c("chrM", "chr2")), filepath = f)
  Rsamtools::indexFa(f)
  FaFile_Circ(Rsamtools::FaFile(f), circularRanges = "chrM")
}
# Ranges on the genome of fasta_circ(), with chrM circular
ranges_circ <- function(x) {
  gr <- GRanges(x)
  seqinfo(gr) <- Seqinfo(seqlevels(gr), seqlengths = rep(10, length(seqlevels(gr))), isCircular = seqlevels(gr) == "chrM")
  gr
}

test_that("FaFile_Circ marks the chromosomes of circularRanges as circular", {
  circ <- isCircular(seqinfo(fasta_circ()))
  expect_true(circ[["chrM"]])
  expect_false(circ[["chr2"]])
})

test_that("getSeq() of a FaFile_Circ reads ranges across the end of a circular chromosome", {
  fa <- fasta_circ()
  expect_equal(as.character(getSeq(fa)), c(chrM = "AAAAAAAAGG", chr2 = "AAAAAAAAGG"))
  grs <- ranges_circ(c("chrM:4-6:+", "chrM:8-11:+", "chrM:8-12:-"))
  expect_no_warning(getSeq(fa, grs[1]))
  expect_equal(as.character(getSeq(fa, grs)), c(chrM = "AAA", chrM = "AGGA", chrM = "TTCCT"))
  # A range across the end of a chromosome that is not circular
  expect_warning(grs <- ranges_circ(c("chrM:4-6:+", "chr2:8-12")), "out-of-bound")
  expect_error(suppressWarnings(getSeq(fa, grs)), "not on circ")
})

test_that("getSeq() of a FaFile_Circ reads ranges that start after the end of a circular chromosome", {
  fa <- fasta_circ()
  # Before, these ranges were dropped, so the sequences did not match the ranges; chrM:11-12 crashed R
  grs <- ranges_circ(c("chrM:19-21:+", "chr2:1-3:+", "chrM:11-12:-", "chrM:28-31:+"))
  expect_equal(as.character(getSeq(fa, grs)), c(chrM = "GGA", chr2 = "AAA", chrM = "TT", chrM = "AGGA"))
})

test_that("extractTranscriptSeqs() of a FaFile_Circ joins exons across the end of a circular chromosome", {
  fa <- fasta_circ()
  txs <- GRangesList(GRanges(c("chrM:4-6:+", "chrM:8-11:+")), GRanges(c("chrM:8-11:-", "chrM:4-6:-")))
  isCircular(seqinfo(txs)) <- TRUE
  seqlengths(txs) <- 10
  expect_equal(extractTranscriptSeqs(fa, txs), DNAStringSet(c("AAAAGGA", "TCCTTTT")))
  # A FaFile cannot
  expect_error(extractTranscriptSeqs(Rsamtools::FaFile(fa$path), txs), "truncated")
})

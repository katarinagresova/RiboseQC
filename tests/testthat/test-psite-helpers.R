# The P-sites of spliced reads. ORFquant uses these functions of RiboseQC too
# (prepare_for_ORFquant()). The P-site is the base at the offset (cutoff) from
# the 5' end of the read, counted on the aligned bases only, so the offset can
# go across an intron. The expected positions are calculated by hand.

# Reads on chr1 with the blocks 100-109 and 160-179 (10M50N20M), or
# 100-104, 115-119 and 130-149 (5M10N5M10N20M)
spliced_reads <- function(strand, cigar = "10M50N20M") {
  GenomicAlignments::GAlignments(seqnames = S4Vectors::Rle("chr1", length(cigar)),
                                 pos = rep(100L, length(cigar)), cigar = cigar,
                                 strand = GenomicRanges::strand(rep(strand, length(cigar))),
                                 seqlengths = c(chr1 = 1000L))
}

test_that("get_ps_fromspliceplus() counts the offset from the start of the read, across the introns", {
  x <- spliced_reads("+")
  expect_identical(start(get_ps_fromspliceplus(x, 5)), 105L)
  expect_identical(start(get_ps_fromspliceplus(x, 9)), 109L)
  expect_identical(start(get_ps_fromspliceplus(x, 10)), 160L)
  expect_identical(start(get_ps_fromspliceplus(x, 12)), 162L)
  expect_identical(start(get_ps_fromspliceplus(spliced_reads("+", "5M10N5M10N20M"), 12)), 132L)
  ps <- get_ps_fromspliceplus(spliced_reads("+", c("10M50N20M", "20M50N10M")), 12)
  expect_identical(start(ps), c(112L, 162L))
  expect_true(all(width(ps) == 1L))
  expect_identical(as.character(strand(ps)), c("+", "+"))
  expect_identical(seqlengths(ps), c(chr1 = 1000L))
})

test_that("get_ps_fromsplicemin() counts the offset from the end of the read, across the introns", {
  x <- spliced_reads("-")
  expect_identical(start(get_ps_fromsplicemin(x, 12)), 167L)
  expect_identical(start(get_ps_fromsplicemin(x, 19)), 160L)
  expect_identical(start(get_ps_fromsplicemin(x, 20)), 109L)
  expect_identical(start(get_ps_fromsplicemin(x, 25)), 104L)
  expect_identical(start(get_ps_fromsplicemin(spliced_reads("-", "5M10N5M10N20M"), 22)), 117L)
  ps <- get_ps_fromsplicemin(spliced_reads("-", c("10M50N20M", "20M50N10M")), 12)
  expect_identical(start(ps), c(117L, 167L))
  expect_true(all(width(ps) == 1L))
  expect_identical(as.character(strand(ps)), c("-", "-"))
  expect_identical(seqlengths(ps), c(chr1 = 1000L))
})

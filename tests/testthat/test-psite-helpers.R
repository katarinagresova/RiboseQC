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

# psites_of_chunk(), for RiboseQC_analysis() and prepare_for_ORFquant(): P_sites_uniq_mm has the
# P-sites of the uniquely mapping reads (mapq more than 50) with mismatches in the MD tag

test_that("psites_of_chunk() finds the reads with mismatches by the MD tag, also in short reads", {
  # 15-nt reads at 101, 201, ..., 601. 1 mismatch in an MD tag of 3 characters (5A9, 9A5) or of 4
  # (0A14), no mismatch, no MD tag, and 1 mismatch in a read that does not map uniquely
  md <- c("5A9", "9A5", "0A14", "15", NA, "5A9")
  x <- GenomicAlignments::GAlignments(seqnames = S4Vectors::Rle("chr1", 6), pos = 100L * (1:6) + 1L,
                                      cigar = rep("15M", 6),
                                      strand = GenomicRanges::strand(c("+", "-", "+", "+", "+", "+")),
                                      seqlengths = c(chr1 = 1000L))
  S4Vectors::mcols(x) <- S4Vectors::DataFrame(mapq = c(255L, 255L, 255L, 255L, 255L, 3L), MD = md)
  junctions <- GenomicRanges::GRanges("chr1", IRanges::IRanges(800, 900), "+", seqlengths = c(chr1 = 1000L))
  res <- psites_of_chunk(x, list(nucl = S4Vectors::DataFrame(read_length = 15, cutoff = 6)), character(),
                         "chr1", c(chr1 = 1000L), junctions)
  # The P-site is 6 nt from the 5' end: 101 + 6 on the plus strand, 215 - 6 on the minus strand
  p_sites <- function(set) sort(start(res[[set]][["15"]]))
  expect_identical(p_sites("P_sites_all"), c(107L, 209L, 307L, 407L, 507L, 607L))
  expect_identical(p_sites("P_sites_uniq"), c(107L, 209L, 307L, 407L, 507L))
  expect_identical(p_sites("P_sites_uniq_mm"), c(107L, 209L, 307L))
})

# unlist(x) gives x back when x is a GRanges or a GAlignments (R/riboseqc.R).
# Without these methods, unlist() of a GRanges stops with an error, and unlist()
# of an empty GRanges gives NULL. RiboseQC_analysis() needs them. ORFquant has no
# copy of the method for GRanges: run_ORFquant() uses this one.

test_that("unlist() of a GRanges or a GAlignments gives it back", {
  gr <- GenomicRanges::GRanges(c("chr1:1-10:+", "chr1:20-30:-"))
  expect_identical(unlist(gr), gr)
  expect_identical(unlist(GenomicRanges::GRanges()), GenomicRanges::GRanges())
  ga <- GenomicAlignments::GAlignments(seqnames = S4Vectors::Rle("chr1", 2L), pos = c(100L, 200L),
                                       cigar = c("10M50N20M", "30M"),
                                       strand = GenomicRanges::strand(c("+", "-")))
  expect_identical(unlist(ga), ga)
})

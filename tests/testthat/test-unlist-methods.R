# RiboseQC 1.4.0 and older had methods of unlist() for GRanges and GAlignments: unlist(x)
# gave x back. S4 methods apply to the whole R session, so these methods also changed
# unlist() in other packages and scripts. Without them, unlist() of a GRanges stops with
# an error, and unlist() of an empty GRanges gives NULL.

test_that("RiboseQC defines no methods of unlist()", {
  expect_length(methods::findMethods("unlist", where = asNamespace("RiboseQC")), 0)
})

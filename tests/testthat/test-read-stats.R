# read_stats of RiboseQC_analysis() on the example data of the package
# (example_analysis() in helper-example.R)

test_that("rld_unq counts in each compartment only its uniquely mapping reads", {
  for (s in c("root", "shoots")) {
    stats <- example_results(s)$read_stats
    # The example genome has the circular chromosomes ChrC and ChrM
    expect_identical(rownames(stats$rld_unq), c("nucl", "ChrC", "ChrM"))
    expect_identical(dimnames(stats$rld_unq), dimnames(stats$rld))
    # The uniquely mapping reads of a compartment are some of its reads
    expect_true(all(stats$rld_unq <= stats$rld), info = s)
  }
})

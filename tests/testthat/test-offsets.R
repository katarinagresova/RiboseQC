# calc_cutoffs_from_profiles() and choose_readlengths() on the 5' profiles of the
# root sample of the example data (example_analysis() in helper-example.R)

test_that("calc_cutoffs_from_profiles() finds the offset of the 28 nt reads, and in-frame signal", {
  profiles <- example_results("root")$profiles_fivepr$five_prime_subcodon
  res <- calc_cutoffs_from_profiles(profiles$nucl[["28"]], length_max = 100)
  expect_identical(res$final_cutoff, 12)
  expect_gt(res$frames_res[["mean_pct_max_frame"]], 80)
  # The same profile with its positions in random order has no frame
  shuffled <- profiles$nucl[["28"]]
  shuffled <- withr::with_seed(1, shuffled[, sample(colnames(shuffled))])
  colnames(shuffled) <- colnames(profiles$nucl[["28"]])
  expect_lt(calc_cutoffs_from_profiles(shuffled, length_max = 100)$frames_res[["mean_pct_max_frame"]], 80)
})

test_that("calc_cutoffs_from_profiles() and choose_readlengths() give the results of the analysis again", {
  # kmeans() in calc_cutoffs_from_profiles() has a fixed seed
  res <- example_results("root")
  profiles <- res$profiles_fivepr$five_prime_subcodon
  expect_identical(calc_cutoffs_from_profiles(profiles$nucl[["28"]], length_max = 28),
                   res$selection_cutoffs$analysis_frame_cutoff$nucl[["28"]])
  expect_identical(choose_readlengths(res$selection_cutoffs$results_cutoffs, choice = "max_coverage",
                                      nt_signals = profiles),
                   res$selection_cutoffs$results_choice)
})

test_that("choose_readlengths() selects a read length with in-frame signal, and none without", {
  profiles <- example_results("root")$profiles_fivepr$five_prime_subcodon
  # One read length (28 nt), with its cutoff and frame preference
  cutoffs <- function(cutoff, frame_preference) {
    list(nucl = matrix(c(cutoff, frame_preference), ncol = 2, dimnames = list("28", c("cutoff", "frame_preference"))))
  }
  res <- calc_cutoffs_from_profiles(profiles$nucl[["28"]], length_max = 100)
  sel <- choose_readlengths(cutoffs(res$final_cutoff, res$frames_res[["mean_pct_max_frame"]]), "max_coverage",
                            profiles["nucl"])
  expect_gt(sel$nucl$data$gain_codons, 0.6)
  expect_equal(nrow(sel$nucl$final_choice), 1)

  shuffled <- profiles$nucl[["28"]]
  shuffled <- withr::with_seed(1, shuffled[, sample(colnames(shuffled))])
  colnames(shuffled) <- colnames(profiles$nucl[["28"]])
  res <- calc_cutoffs_from_profiles(shuffled, length_max = 100)
  sel <- choose_readlengths(cutoffs(res$final_cutoff, res$frames_res[["mean_pct_max_frame"]]), "max_coverage",
                            profiles["nucl"])
  expect_lt(sel$nucl$data$gain_codons, 0)
  expect_equal(nrow(sel$nucl$final_choice), 0)
})

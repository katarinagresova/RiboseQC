# calc_cutoffs_from_profiles() and choose_readlengths() on the 5' profiles of ChrC of the
# shoots sample of the example data (example_analysis() in helper-example.R). In nucl,
# no read length of the example has an offset. The tests at the end use hand-made 5'
# profiles

test_that("calc_cutoffs_from_profiles() finds the offset of the 29 nt reads, and in-frame signal", {
  profiles <- example_results("shoots")$profiles_fivepr$five_prime_subcodon
  res <- calc_cutoffs_from_profiles(profiles$ChrC[["29"]], length_max = 100)
  expect_identical(res$final_cutoff, 13)
  expect_gt(res$frames_res[["mean_pct_max_frame"]], 70)
  # The same profile with its positions in random order has no frame
  shuffled <- profiles$ChrC[["29"]]
  shuffled <- withr::with_seed(1, shuffled[, sample(colnames(shuffled))])
  colnames(shuffled) <- colnames(profiles$ChrC[["29"]])
  expect_lt(calc_cutoffs_from_profiles(shuffled, length_max = 100)$frames_res[["mean_pct_max_frame"]], 70)
})

test_that("calc_cutoffs_from_profiles() and choose_readlengths() give the results of the analysis again", {
  # kmeans() in calc_cutoffs_from_profiles() has a fixed seed
  res <- example_results("shoots")
  profiles <- res$profiles_fivepr$five_prime_subcodon
  expect_identical(calc_cutoffs_from_profiles(profiles$ChrC[["29"]], length_max = 29),
                   res$selection_cutoffs$analysis_frame_cutoff$ChrC[["29"]])
  expect_identical(choose_readlengths(res$selection_cutoffs$results_cutoffs, choice = "max_coverage",
                                      nt_signals = profiles),
                   res$selection_cutoffs$results_choice)
})

test_that("choose_readlengths() selects a read length with in-frame signal, and none without", {
  profiles <- example_results("shoots")$profiles_fivepr$five_prime_subcodon
  # One read length (29 nt), with its cutoff and frame preference
  cutoffs <- function(cutoff, frame_preference) {
    list(ChrC = matrix(c(cutoff, frame_preference), ncol = 2, dimnames = list("29", c("cutoff", "frame_preference"))))
  }
  res <- calc_cutoffs_from_profiles(profiles$ChrC[["29"]], length_max = 100)
  sel <- choose_readlengths(cutoffs(res$final_cutoff, res$frames_res[["mean_pct_max_frame"]]), "max_coverage",
                            profiles["ChrC"])
  expect_gt(sel$ChrC$data$gain_codons, 0.3)
  expect_equal(nrow(sel$ChrC$final_choice), 1)

  shuffled <- profiles$ChrC[["29"]]
  shuffled <- withr::with_seed(1, shuffled[, sample(colnames(shuffled))])
  colnames(shuffled) <- colnames(profiles$ChrC[["29"]])
  res <- calc_cutoffs_from_profiles(shuffled, length_max = 100)
  sel <- choose_readlengths(cutoffs(res$final_cutoff, res$frames_res[["mean_pct_max_frame"]]), "max_coverage",
                            profiles["ChrC"])
  expect_lt(sel$ChrC$data$gain_codons, 0)
  expect_equal(nrow(sel$ChrC$final_choice), 0)
})

# Hand-made 5' profiles, where the offset and the frame are known. The columns are those of
# five_prime_subcodon: 5_UTR_1..50, CDS_1..99 and 3_UTR_1..50. calc_cutoffs_from_profiles()
# and choose_readlengths() take them as one sequence: CDS_1 is the first nt of the start
# codon, and 5_UTR_50 is the nt before it.
# A profile of n_tx transcripts, for reads with the offset `offset` (one value, or one for
# each transcript): on each transcript, `start` reads with the P-site on the start codon,
# and `reads` reads on each codon of `codons` (codon 1 is the codon after the start codon).
# With frame = FALSE, the 5' ends of the reads of a codon are on each of its 3 nt
fivepr_profile <- function(offset, n_tx = 20, start = 10, codons = 1:20, reads = 1, frame = TRUE) {
  cols <- c(paste0("5_UTR_", 1:50), paste0("CDS_", 1:99), paste0("3_UTR_", 1:50))
  prof <- matrix(0, n_tx, length(cols), dimnames = list(paste0("tx", seq_len(n_tx)), cols))
  offset <- rep_len(offset, n_tx)
  for (t in seq_len(n_tx)) {
    # the column of the 5' end of a read with the P-site on the first nt of the start codon
    atg <- which(cols == "CDS_1") - offset[t]
    ends <- atg + 3 * codons
    if (!frame) ends <- c(ends, ends + 1, ends + 2)
    prof[t, atg] <- start
    prof[t, ends] <- reads
  }
  S4Vectors::DataFrame(prof, check.names = FALSE)
}

# The cutoff and the frame preference of each read length of `profiles`, the profiles of one
# compartment, as RiboseQC_analysis() gives them to choose_readlengths()
cutoffs_of <- function(profiles) {
  t(sapply(profiles, function(p) {
    res <- calc_cutoffs_from_profiles(p, length_max = 100)
    c(cutoff = res$final_cutoff, frame_preference = res$frames_res[["mean_pct_max_frame"]])
  }))
}

test_that("calc_cutoffs_from_profiles() finds the offset of hand-made profiles, in each frame", {
  for (offset in c(0, 12, 13, 14, 24)) {
    res <- calc_cutoffs_from_profiles(fivepr_profile(offset), length_max = 100)
    expect_identical(res$final_cutoff, offset)
    expect_equal(res$frames_res, c(max_frame_in_pctORF = 100, mean_pct_max_frame = 100))
  }
})

test_that("calc_cutoffs_from_profiles() gives the most frequent offset, with each of the 8 methods", {
  # 30 transcripts with the offset 12, and 20 with 4 other offsets in the same frame. The
  # methods with kmeans() need more than 4 different values
  prof <- fivepr_profile(rep(c(12, 15, 9, 18, 6), c(30, 10, 6, 3, 1)), n_tx = 50)
  res <- calc_cutoffs_from_profiles(prof, length_max = 100)
  expect_identical(res$final_cutoff, 12)
  expect_identical(unname(res$cutoffs), rep(12, 8))
})

test_that("calc_cutoffs_from_profiles() finds the offset from 24 nt before the start codon to its first nt, and up to length_max", {
  res <- calc_cutoffs_from_profiles(fivepr_profile(30), length_max = 100)
  expect_identical(colnames(res$profiles_start), as.character(-24:0))
  expect_lte(res$final_cutoff, 24)
  res <- calc_cutoffs_from_profiles(fivepr_profile(15), length_max = 12)
  expect_identical(colnames(res$profiles_start), as.character(-12:0))
  expect_lte(res$final_cutoff, 12)
})

test_that("calc_cutoffs_from_profiles() gives NULL if fewer than 5 transcripts have 5' ends in the CDS, or near the start codon", {
  expect_null(calc_cutoffs_from_profiles(fivepr_profile(12, n_tx = 4), length_max = 100))
  expect_identical(calc_cutoffs_from_profiles(fivepr_profile(12, n_tx = 5), length_max = 100)$final_cutoff, 12)
  # The transcripts without reads do not count
  prof <- rbind(fivepr_profile(12, n_tx = 4), fivepr_profile(12, n_tx = 16, start = 0, codons = integer()))
  expect_null(calc_cutoffs_from_profiles(prof, length_max = 100))
  # 20 transcripts with 5' ends before the start codon, but only 4 with 5' ends in the CDS
  prof <- rbind(fivepr_profile(12, n_tx = 4), fivepr_profile(12, n_tx = 16, codons = integer()))
  expect_null(calc_cutoffs_from_profiles(prof, length_max = 100))
  # 20 transcripts with 5' ends in the CDS, but only 4 with 5' ends from 24 nt before the
  # start codon to its first nt
  prof <- rbind(fivepr_profile(12, n_tx = 4), fivepr_profile(12, n_tx = 16, start = 0, codons = 10:20))
  expect_null(calc_cutoffs_from_profiles(prof, length_max = 100))
})

test_that("calc_cutoffs_from_profiles() gives the frame preference of the 5' ends in the CDS", {
  res <- calc_cutoffs_from_profiles(fivepr_profile(12, frame = FALSE), length_max = 100)
  expect_identical(res$final_cutoff, 12)
  expect_equal(res$frames_res[["mean_pct_max_frame"]], 100 / 3)
  expect_equal(unname(res$frames), matrix(1 / 3, 20, 3))
})

test_that("choose_readlengths() selects the read lengths with in-frame signal, with each method", {
  profiles <- list(nucl = list("30" = fivepr_profile(12, frame = FALSE), "28" = fivepr_profile(12),
                               "29" = fivepr_profile(13), all = fivepr_profile(12)))
  cutoffs <- list(nucl = cutoffs_of(profiles$nucl))
  for (choice in c("max_coverage", "max_inframe")) {
    sel <- choose_readlengths(cutoffs, choice, profiles)$nucl
    expect_identical(sel$final_choice$read_length, c("28", "29"))
    expect_identical(sel$final_choice$cutoff, c(12, 13))
  }
  sel <- choose_readlengths(cutoffs, "all", profiles)$nucl
  expect_identical(sel$final_choice$read_length, c("28", "29", "30"))
  # The values of each read length but "all", from the highest frame preference to the lowest
  expect_identical(sel$data$read_length, c("28", "29", "30"))
  expect_equal(sel$data$gain_codons[1:2], c(1, 1))
  expect_lt(sel$data$gain_codons[3], 0)
  expect_identical(sel$data$max_coverage, c(TRUE, TRUE, FALSE))
})

test_that("choose_readlengths(choice = \"max_inframe\") selects only the read lengths that add in-frame codons", {
  # The in-frame reads of 29 nt are on the codons where 28 nt has more than 1 read, so
  # they are masked. 29 nt has also reads on 5 other codons, with the offset 12: with its
  # cutoff 13, they are out of frame
  in_frame <- as.matrix(fivepr_profile(13))
  out_of_frame <- as.matrix(fivepr_profile(12, start = 0, codons = 21:25))
  profiles <- list(nucl = list("28" = fivepr_profile(12, reads = 2),
                               "29" = S4Vectors::DataFrame(in_frame + out_of_frame, check.names = FALSE)))
  cutoffs <- list(nucl = cutoffs_of(profiles$nucl))
  expect_identical(cutoffs$nucl[, "cutoff"], c("28" = 12, "29" = 13))
  sel <- choose_readlengths(cutoffs, "max_coverage", profiles)$nucl
  expect_identical(sel$final_choice$read_length, c("28", "29"))
  sel <- choose_readlengths(cutoffs, "max_inframe", profiles)$nucl
  expect_identical(sel$final_choice$read_length, "28")
  expect_gt(sel$data$gain_codons[2], 0)
  expect_lt(sel$data$gain_new_codons[2], 0)
})

test_that("choose_readlengths() stops for an unknown method, and for the profiles of other compartments", {
  profiles <- list(nucl = list("28" = fivepr_profile(12)))
  cutoffs <- list(nucl = cutoffs_of(profiles$nucl))
  expect_error(choose_readlengths(cutoffs, "best", profiles), "choice must be one of")
  expect_error(choose_readlengths(list(ChrC = cutoffs$nucl), "max_coverage", profiles),
               "names(summary_data) == names(nt_signals)", fixed = TRUE)
})

test_that("choose_readlengths() gives no result for a compartment with only the read length \"all\"", {
  profiles <- list(nucl = list(all = fivepr_profile(12)), ChrC = list("28" = fivepr_profile(12)))
  cutoffs <- lapply(profiles, cutoffs_of)
  expect_identical(names(choose_readlengths(cutoffs, "max_coverage", profiles)), "ChrC")
})

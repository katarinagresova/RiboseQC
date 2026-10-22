# The functions that take the data of the report from the results of RiboseQC_analysis(), on the
# bundled results of the root BAM file (inst/ext_data/root_results_RiboseQC)

root_results <- function() get(load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC")))

test_that("get_metagene_data() and get_codon_usage_data() name comp in the error for a wrong comp", {
  res_all <- root_results()
  # Before, the error said "Variable res"
  expect_error(get_metagene_data(res_all, "profiles_fivepr", "five_prime_subcodon", "no_comp"),
               "Variable comp takes only: .*nucl")
  expect_error(get_codon_usage_data(res_all, "Codon_counts", "no_comp", "28"), "Variable comp takes only: .*nucl")
})

# The functions that take the data of the report from the results of RiboseQC_analysis(), on the
# bundled results of the root and shoots BAM files (inst/ext_data/root_results_RiboseQC and
# shoots_results_RiboseQC). The root sample has no P-sites, so no codon tables

root_results <- function() get(load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC")))
shoots_results <- function() get(load(system.file("ext_data", "shoots_results_RiboseQC", package = "RiboseQC")))

test_that("get_metagene_data() and get_codon_usage_data() name comp in the error for a wrong comp", {
  res_all <- root_results()
  # Before, the error said "Variable res"
  expect_error(get_metagene_data(res_all, "profiles_fivepr", "five_prime_subcodon", "no_comp"),
               "Variable comp takes only: .*nucl")
  expect_error(get_codon_usage_data(shoots_results(), "Codon_counts", "no_comp", "28"), "Variable comp takes only: ChrC")
})

test_that("get_top50_cds_genes() and get_top50_all_genes() give at most 50 and 100 rows, without rows of NA", {
  res_all <- root_results()
  # The results with only the first n genes
  top <- function(n) {
    r <- res_all
    r$read_stats$counts_cds_genes <- r$read_stats$counts_cds_genes[seq_len(n), ]
    r$read_stats$counts_all_genes <- r$read_stats$counts_all_genes[seq_len(n), ]
    list(cds = get_top50_cds_genes(list(root = r))$root, all = get_top50_all_genes(list(root = r))$root)
  }
  all_genes <- top(nrow(res_all$read_stats$counts_cds_genes))
  expect_equal(c(nrow(all_genes$cds), nrow(all_genes$all)), c(50, 100))
  expect_identical(all_genes$cds$reads, head(sort(res_all$read_stats$counts_cds_genes$reads, decreasing = TRUE), 50))
  # Before, fewer genes gave rows of NA, up to 50 and 100 rows
  few <- top(10)
  expect_equal(c(nrow(few$cds), nrow(few$all)), c(10, 10))
  expect_false(anyNA(few$cds$reads))
  expect_false(anyNA(few$all$reads))
})

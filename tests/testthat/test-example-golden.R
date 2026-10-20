# The results of RiboseQC_analysis() on the example data of the package
# (example_analysis() in helper-example.R) must be those of golden/, which are
# what RiboseQC gave when the files were made. For each sample: the offsets
# (P_sites_offsets.tsv), the statistics of the selected read lengths
# (P_sites_calcs), the reads of each read length and of each region in each
# compartment (read_stats$rld and read_stats$positions), the reads, RPKM and TPM
# of each gene (read_stats$counts_all_genes), and the P-sites
# (P_sites_plus.bedgraph and P_sites_minus.bedgraph).
#
# If a change that must not change results makes this test fail, the change is
# wrong. If the results must change, make the golden files again and review
# their diff, from the package directory:
#   RIBOSEQC_UPDATE_GOLDEN=1 Rscript -e 'testthat::test_local()'

golden_names <- c("P_sites_offsets.tsv", "P_sites_calcs", "read_lengths.tsv", "regions.tsv", "genes.tsv",
                  "P_sites_plus.bedgraph", "P_sites_minus.bedgraph")

# The results of `sample` as files with the names of golden_names; their paths
result_files <- function(sample) {
  prefix <- example_analysis()$prefix[[sample]]
  stats <- example_results(sample)$read_stats
  files <- stats::setNames(paste0(prefix, "_", golden_names), golden_names)
  write_tsv <- function(df, f) utils::write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
  write_tsv(data.frame(compartment = rownames(stats$rld), stats$rld, check.names = FALSE), files[["read_lengths.tsv"]])
  write_tsv(data.frame(region = rownames(stats$positions), stats$positions, check.names = FALSE), files[["regions.tsv"]])
  write_tsv(data.frame(gene_id = rownames(stats$counts_all_genes), as.data.frame(stats$counts_all_genes),
                       check.names = FALSE), files[["genes.tsv"]])
  files
}
# The golden files of `sample`; the bedgraph files are compressed
golden_files <- function(sample) {
  f <- test_path("golden", paste0(sample, "_", golden_names))
  f <- ifelse(endsWith(f, ".bedgraph"), paste0(f, ".gz"), f)
  stats::setNames(f, golden_names)
}

# A result file as a data.frame; the bedgraph files have no header
read_result <- function(f) {
  utils::read.delim(f, header = !grepl("\\.bedgraph", f), quote = "", stringsAsFactors = FALSE, check.names = FALSE)
}

for (s in c("root", "shoots")) {
  test_that(paste("the results of the", s, "sample are those of the golden files"), {
    new <- result_files(s)
    old <- golden_files(s)
    if (nzchar(Sys.getenv("RIBOSEQC_UPDATE_GOLDEN"))) {
      dir.create(test_path("golden"), showWarnings = FALSE)
      for (n in golden_names) {
        if (endsWith(old[[n]], ".gz")) {
          gz <- gzfile(old[[n]], "w")
          writeLines(readLines(new[[n]]), gz)
          close(gz)
        } else {
          file.copy(new[[n]], old[[n]], overwrite = TRUE)
        }
      }
    }
    # The numbers to 6 digits
    for (n in golden_names) {
      expect_equal(read_result(new[[n]]), read_result(old[[n]]), tolerance = 1e-6, info = n)
    }
    # Golden files made from an empty result would pass the comparison
    expect_gt(nrow(read_result(new[["P_sites_offsets.tsv"]])), 0)
    expect_gt(nrow(read_result(new[["P_sites_plus.bedgraph"]])), 0)
  })
}

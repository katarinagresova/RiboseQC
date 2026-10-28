# The results of RiboseQC_analysis() on the example data of the package
# (example_analysis() in helper-example.R) must be those of golden/, which are
# what RiboseQC gave when the files were made. For each sample: the offsets
# (P_sites_offsets.tsv), the statistics of the selected read lengths
# (P_sites_calcs), the reads of each read length and of each region in each
# compartment (read_stats$rld and read_stats$positions, in read_lengths.tsv and
# regions.tsv), the reads, RPKM and TPM of each gene
# (read_stats$counts_all_genes, in genes.tsv), and the P-sites
# (P_sites_plus.bedgraph and P_sites_minus.bedgraph). These are the files of
# RiboseQC_analysis(). No read length of root has an offset, so root has no
# P-sites and no P-site files.
#
# If a change that must not change results makes this test fail, the change is
# wrong. If the results must change, make the golden files again and review
# their diff, from the package directory:
#   RIBOSEQC_UPDATE_GOLDEN=1 Rscript -e 'testthat::test_local()'

golden_names <- c("P_sites_offsets.tsv", "P_sites_calcs", "read_lengths.tsv", "regions.tsv", "genes.tsv",
                  "P_sites_plus.bedgraph", "P_sites_minus.bedgraph")

# The paths of the result files of `sample` with the names of golden_names
result_files <- function(sample) {
  stats::setNames(paste0(example_analysis()$prefix[[sample]], "_", golden_names), golden_names)
}
# The golden files of `sample`; the bedgraph files are compressed
golden_files <- function(sample) {
  f <- test_path("golden", paste0(sample, "_", golden_names))
  f <- ifelse(endsWith(f, ".bedgraph"), paste0(f, ".gz"), f)
  stats::setNames(f, golden_names)
}

# A result file as a data.frame; the bedgraph files have no header. An empty file
# (P_sites_calcs without offsets) gives a data.frame without columns
read_result <- function(f) {
  if (!any(nzchar(readLines(f)))) return(data.frame())
  utils::read.delim(f, header = !grepl("\\.bedgraph", f), quote = "", stringsAsFactors = FALSE, check.names = FALSE)
}

for (s in c("root", "shoots")) {
  test_that(paste("the results of the", s, "sample are those of the golden files"), {
    new <- result_files(s)
    old <- golden_files(s)
    files <- if (s == "root") grep("bedgraph", golden_names, value = TRUE, invert = TRUE) else golden_names
    if (nzchar(Sys.getenv("RIBOSEQC_UPDATE_GOLDEN"))) {
      dir.create(test_path("golden"), showWarnings = FALSE)
      for (n in files) {
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
    for (n in files) {
      expect_equal(read_result(new[[n]]), read_result(old[[n]]), tolerance = 1e-6, info = n)
    }
    if (s == "root") {
      expect_identical(nrow(read_result(new[["P_sites_offsets.tsv"]])), 0L)
      expect_false(file.exists(new[["P_sites_plus.bedgraph"]]))
    } else {
      # Golden files made from an empty result would pass the comparison
      expect_gt(nrow(read_result(new[["P_sites_offsets.tsv"]])), 0)
      expect_gt(nrow(read_result(new[["P_sites_plus.bedgraph"]])), 0)
    }
  })
}

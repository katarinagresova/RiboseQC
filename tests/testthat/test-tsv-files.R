# The tab-separated files of RiboseQC_analysis(write_TSV_files = TRUE), the default, on
# the example data (example_analysis() in helper-example.R): tables of res_all, and a
# summary of each compartment. test-example-golden.R compares read_lengths.tsv,
# regions.tsv and genes.tsv with golden files

# A tab-separated file as a data.frame; the columns `text` as they are written, also "NA".
# An empty file (P_sites_calcs without offsets) gives a data.frame without columns
read_tsv <- function(f, text = character()) {
  if (!any(nzchar(readLines(f)))) return(data.frame())
  df <- utils::read.delim(f, quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  if (length(text) > 0) {
    df[text] <- utils::read.delim(f, quote = "", colClasses = "character", na.strings = character(),
                                  check.names = FALSE)[text]
  }
  df
}

# The summary of `res_all` with the offsets of P_sites_offsets.tsv `offs` and the table of
# P_sites_calcs `calcs`: for each compartment, the percentages of its reads on the CDS and
# with a read length that has an offset, and the frame preference of these read lengths,
# weighted by their reads
expected_summary <- function(res_all, offs, calcs) {
  rld <- res_all$read_stats$rld
  do.call(rbind, lapply(rownames(rld), function(comp) {
    reads <- sum(rld[comp, ])
    rls <- offs$read_length[offs$compartment == comp]
    rl_reads <- vapply(paste0("reads_", rls), function(n) if (n %in% names(rld)) rld[comp, n] else 0, 0)
    fp <- calcs$frame_preference[match(paste(comp, rls), paste(calcs$comp, calcs$read_length))]
    has_fp <- !is.na(fp) & rl_reads > 0
    data.frame(compartment = comp, reads = reads,
               pct_cds = round(res_all$read_stats$positions["cds", paste0("reads_", comp)] / reads * 100, 4),
               pct_P_sites = round(sum(rl_reads) / reads * 100, 4),
               frame_preference = if (any(has_fp)) round(stats::weighted.mean(fp[has_fp], rl_reads[has_fp]), 4) else NA_real_)
  }))
}

for (s in c("root", "shoots")) {
  test_that(paste("RiboseQC_analysis() writes the tables of the", s, "sample as tab-separated files"), {
    prefix <- example_analysis()$prefix[[s]]
    res_all <- example_results(s)

    expect_identical(readLines(paste0(prefix, "_P_sites_calcs.tsv")), readLines(paste0(prefix, "_P_sites_calcs")))

    top <- read_tsv(paste0(prefix, "_top50_positions.tsv"), text = c("seq", "region", "gene", "gene_biotype"))
    pos <- res_all$sequence_analysis
    expect_identical(names(top), c("seqnames", "start", "end", "width", "strand", names(mcols(pos))))
    expect_identical(nrow(top), 50L)
    expect_identical(top$seqnames, as.character(seqnames(pos)))
    expect_identical(top$start, start(pos))
    expect_identical(top$end, end(pos))
    expect_identical(top$strand, as.character(strand(pos)))
    expect_equal(top$score, pos$score)
    expect_equal(top$len_adj, pos$len_adj)
    expect_equal(top$pct, pos$pct)
    expect_identical(top$seq, unname(as.character(pos$seq)))
    for (n in c("region", "gene", "gene_biotype")) {
      expect_identical(top[[n]], vapply(as.list(mcols(pos)[[n]]), paste, "", collapse = ","), info = n)
    }

    summ <- read_tsv(paste0(prefix, "_summary.tsv"))
    # all NA for root: read.delim() reads the column as logical
    summ$frame_preference <- as.numeric(summ$frame_preference)
    # example_analysis() gives no sample_names: the default names
    expect_identical(summ$sample, rep(c(root = "sample1", shoots = "sample2")[[s]], 3))
    expected <- expected_summary(res_all, read_tsv(paste0(prefix, "_P_sites_offsets.tsv")),
                                 read_tsv(paste0(prefix, "_P_sites_calcs")))
    expect_equal(summ[-1], expected)
    # root has no offsets; shoots has offsets only for ChrC
    expect_identical(!is.na(summ$frame_preference), c(FALSE, s == "shoots", FALSE))
  })
}

test_that("the summary has the reads of all read lengths with an offset, and the frame preference of those with one", {
  res_all <- example_results("shoots")
  finn <- res_all$summary_P_sites
  prefix <- file.path(withr::local_tempdir(), "shoots")
  summary_of <- function(...) {
    write_tsv_files(res_all, ..., dest_name = prefix, sample_name = "shoots")
    read_tsv(paste0(prefix, "_summary.tsv"))
  }
  rld <- res_all$read_stats$rld

  # 33 nt has an offset but no frame preference, as with rescue_all_rls or offsets_df; 99 nt has no reads
  summ <- summary_of(finn, data.frame(read_length = c(28, 33, 99), cutoff = 12, compartment = "ChrC"))
  expect_equal(summ$pct_P_sites[2], round((rld["ChrC", "reads_28"] + rld["ChrC", "reads_33"]) / sum(rld["ChrC", ]) * 100, 4))
  expect_equal(summ$frame_preference[2], round(finn$frame_preference[finn$read_length == "28" & finn$comp == "ChrC"], 4))

  # no read length has an offset
  summ <- summary_of(NULL, data.frame(read_length = numeric(), cutoff = numeric(), compartment = character()))
  expect_equal(summ$pct_P_sites, c(0, 0, 0))
  expect_true(all(is.na(summ$frame_preference)))
  expect_identical(readLines(paste0(prefix, "_P_sites_calcs.tsv")), "")

  # a compartment without reads
  res_all$read_stats$rld["ChrM", ] <- 0
  res_all$read_stats$positions[, "reads_ChrM"] <- 0
  summ <- summary_of(finn, data.frame(read_length = 28, cutoff = 12, compartment = "ChrC"))
  expect_identical(summ$reads[3], 0L)
  expect_true(is.na(summ$pct_cds[3]) && is.na(summ$pct_P_sites[3]) && is.na(summ$frame_preference[3]))
})

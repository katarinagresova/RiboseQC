# The 5' profiles of the transcripts of RiboseQC_analysis() on the example data
# (example_analysis() in helper-example.R)

test_that("each row of the 5' profiles has the CDS windows of the transcript of the row", {
  # Before, the 5'UTR and 3'UTR windows were in the order of the transcript names and the
  # CDS windows in another order, so some rows had the CDS windows of another transcript
  res <- example_results("root")
  nucl <- res$profiles_fivepr$five_prime_subcodon$nucl
  prof <- as.matrix(nucl[["all"]])
  txs <- rownames(prof)
  ends <- unlist(res$read_stats$reads_pos1[setdiff(names(nucl), "all")], use.names = FALSE)
  cds <- get(load(example_analysis()$annotation))$cds_txs[txs]
  mp <- GenomicFeatures::mapToTranscripts(ends, cds)
  seqlengths(mp) <- sum(width(cds))[seqlevels(mp)]
  cov <- coverage(mp, weight = ends$score[mp$xHits])
  # the 5' ends on the first 33 nt of the CDS of each transcript
  first33 <- t(vapply(txs, function(tx) as.vector(cov[[tx]])[1:33], numeric(33)))
  expect_gt(sum(first33), 0)
  expect_equal(unname(prof[, paste0("CDS_", 1:33)]), unname(first33))
})

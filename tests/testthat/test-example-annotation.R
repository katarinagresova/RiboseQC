# The annotation of the example data of the package (example_analysis() in
# helper-example.R): 309 genes on the 5 nuclear chromosomes, the chloroplast genome
# (ChrC) and the mitochondrial genome (ChrM)

test_that("prepare_annotation_files() makes the annotation of the example data", {
  annotation <- example_analysis()$annotation
  # The TxDb, and the tables of export_bed_tables_TxDb = TRUE
  expect_true(all(file.exists(file.path(dirname(annotation), c(basename(annotation), "simp_arab.gtf.gz_TxDb",
                                                               "cds_txs_coords_similbed.bed", "table_gene_tx_IDs")))))
  env <- new.env()
  suppressMessages(load_annotation(annotation, envir = env))
  ann <- env$GTF_annotation
  expect_length(ann$genes, 309)
  expect_length(ann$txs, 410)
  expect_equal(c(table(as.character(seqnames(ann$genes)))[c("ChrC", "ChrM")]), c(ChrC = 88, ChrM = 122))
  # The default circ_chroms make ChrC and ChrM circular. All chromosomes have the standard genetic code
  expect_setequal(names(which(isCircular(ann$seqinfo))), c("ChrC", "ChrM"))
  expect_true(all(ann$genetic_codes$genetic_code == "1"))
  expect_s4_class(env$genome_seq, "FaFile_Circ")
})

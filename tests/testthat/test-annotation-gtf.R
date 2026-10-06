# The annotation of GTF files that are not GENCODE's. The GTF of the example
# (GENCODE's attribute names, all biotypes, no gene on two strands) is changed in
# the ways that other GTFs differ from it. The annotation must be that of the
# example, or the one that NEWS.md describes. From ORFquant's tests
# (tests/testthat/test-annotation-gtf.R), without the test of run_ORFquant().

work <- tempfile("riboseqc_gtf_")
dir.create(work)
gtf_example <- readLines(gzfile(example_file("chr22_example.gtf.gz")))
fasta <- example_file("chr22_example.fa.gz")

# GenomicFeatures reads a start or stop codon as a CDS too
coding_types <- c("CDS", "start_codon", "stop_codon")
line_type <- function(lines) sub("^[^\t]*\t[^\t]*\t([^\t]*)\t.*", "\\1", lines)
line_attr <- function(lines, key) {
  has <- grepl(paste0(key, ' "'), lines, fixed = TRUE)
  ifelse(has, sub(paste0('.*', key, ' "([^"]*)".*'), "\\1", lines), NA_character_)
}
# The lines without the attribute `key`, on the `lines` selected by `where`
drop_attr <- function(lines, key, where = rep(TRUE, length(lines))) {
  lines[where] <- sub(paste0(key, ' "[^"]*"; ?'), "", lines[where])
  lines
}

# Writes `lines` as a GTF file, prepares its annotation, and gives the path of the
# annotation file; `annotate()` gives the annotation
annotation_of <- function(lines, name) {
  dir <- file.path(work, name)
  dir.create(dir)
  gtf <- file.path(dir, "genes.gtf")
  writeLines(lines, gtf)
  suppressMessages(suppressWarnings(
    prepare_annotation_files(annotation_directory = dir, gtf_file = gtf,
                             genome_seq = fasta, forge_BSgenome = FALSE)))
  file.path(dir, "genes.gtf_Rannot")
}
annotate <- function(lines, name) get(load(annotation_of(lines, name)))

ann_genc <- annotate(gtf_example, "genc")
tx_ids <- unique(stats::na.omit(line_attr(gtf_example, "transcript_id")))

test_that("the annotation of the example has one row for each transcript", {
  expect_equal(ann_genc$trann$transcript_id, tx_ids)
  expect_true(all(ann_genc$trann$transcript_biotype %in%
                    c("protein_coding", "protein_coding_CDS_not_defined", "retained_intron",
                      "nonsense_mediated_decay", "non_stop_decay")))
})

# Names of attributes, as in the GTFs of Ensembl, NCBI, UCSC, gffread and StringTie
test_that("ids, names and biotypes are read by the name of their attribute", {
  # other names of the attributes, also mixed with GENCODE's: those of Ensembl, NCBI and StringTie GTFs
  renames <- list(c(gene_type = "gene_biotype", gene_name = "gene_symbol", transcript_type = "transcript_biotype"),
                  c(gene_name = "gene", transcript_type = "transcript_biotype"),
                  c(gene_name = "ref_gene_name"))
  for (i in seq_along(renames)) {
    gtf <- gtf_example
    for (old in names(renames[[i]])) gtf <- sub(paste0(old, ' "'), paste0(renames[[i]][[old]], ' "'), gtf)
    expect_equal(annotate(gtf, paste0("renamed", i))$trann, ann_genc$trann, info = paste(renames[[i]], collapse = ", "))
  }
})

test_that("a biotype or a name that only some lines of a gene or transcript have is found", {
  # as in the GTF of gffread, which writes them on the transcript lines only
  other <- !line_type(gtf_example) %in% c("gene", "transcript")
  gtf <- gtf_example
  for (key in c("gene_type", "gene_name", "transcript_type")) gtf <- drop_attr(gtf, key, other)
  expect_equal(annotate(gtf, "tx_lines")$trann, ann_genc$trann)
})

test_that("a GTF without gene names has the gene name no_name", {
  trann <- annotate(drop_attr(gtf_example, "gene_name"), "no_names")$trann
  expect_equal(trann$gene_name, rep("no_name", nrow(trann)))
  expect_equal(trann$gene_biotype, ann_genc$trann$gene_biotype)
})

test_that("a GTF without CDS lines is refused", {
  msg <- tryCatch(annotation_of(gtf_example[!line_type(gtf_example) %in% coding_types], "no_cds"),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "no CDS lines")
})

test_that("a transcript with a CDS shorter than 3 nt is not a coding transcript", {
  # a copy of a coding transcript, with one CDS line of 2 nt instead of its CDS
  cds_lines <- gtf_example[line_type(gtf_example) == "CDS"]
  tx <- line_attr(cds_lines, "transcript_id")[1]
  copy <- gtf_example[line_attr(gtf_example, "transcript_id") %in% tx &
                        !line_type(gtf_example) %in% c(coding_types, "UTR")]
  cds <- strsplit(cds_lines[line_attr(cds_lines, "transcript_id") == tx][1], "\t")[[1]]
  cds[5] <- as.integer(cds[4]) + 1L
  copy <- gsub(tx, paste0(tx, "_short"), c(copy, paste(cds, collapse = "\t")), fixed = TRUE)
  ann <- annotate(c(gtf_example, copy), "short_cds")
  expect_true(paste0(tx, "_short") %in% ann$trann$transcript_id)
  expect_true(paste0(tx, "_short") %in% names(ann$exons_txs))
  # the coding transcripts and their CDS are those of the example
  cds_ranges <- function(a) { x <- unlist(a$cds_txs); mcols(x) <- NULL; x }
  expect_equal(cds_ranges(ann), cds_ranges(ann_genc))
  expect_equal(ann$cds_txs_coords$reprentative_mostcommon, ann_genc$cds_txs_coords$reprentative_mostcommon)
})

# Biotypes
test_that("the biotype mRNA of NCBI is protein_coding", {
  gtf <- sub('transcript_type "protein_coding"', 'transcript_type "mRNA"', gtf_example)
  trann <- annotate(gtf, "mrna")$trann
  expect_false("mRNA" %in% trann$transcript_biotype)
  expect_equal(trann, ann_genc$trann)
})

test_that("without biotypes, a transcript with a CDS and its gene are protein_coding", {
  # RANBP1 has no CDS lines here (and no start and stop codons): it is not coding
  ranbp1 <- "ENSG00000099901.18"
  gtf <- drop_attr(drop_attr(gtf_example, "gene_type"), "transcript_type")
  gtf <- gtf[!(line_type(gtf) %in% coding_types & line_attr(gtf, "gene_id") == ranbp1)]
  trann <- annotate(gtf, "no_types")$trann
  has_cds <- trann$transcript_id %in% line_attr(gtf, "transcript_id")[line_type(gtf) == "CDS"]
  expect_equal(trann$transcript_biotype, ifelse(has_cds, "protein_coding", "no_type"))
  expect_equal(trann$gene_biotype, ifelse(trann$gene_id == ranbp1, "no_type", "protein_coding"))
  # also those that GENCODE calls nonsense_mediated_decay, which have a CDS
  nmd <- ann_genc$trann$transcript_biotype == "nonsense_mediated_decay" & ann_genc$trann$gene_id != ranbp1
  expect_true(any(nmd))
  expect_true(all(trann$transcript_biotype[nmd] == "protein_coding"))
})

# Genes that genes() drops, as the PAR genes of UCSC's GTFs, which are on chrX and chrY.
# DDT with a copy of its transcript ENST00000350608.7 on the other strand: the gene
# is on both strands, so genes() drops it
ddt <- "ENSG00000099977.16"
ddt_tx <- "ENST00000350608.7"
gtf_two_strands <- local({
  copy <- gtf_example[line_attr(gtf_example, "transcript_id") %in% ddt_tx]
  copy <- sub("\t-\t", "\t+\t", sub(ddt_tx, paste0(ddt_tx, "_rev"), copy, fixed = TRUE))
  c(gtf_example, copy)
})

test_that("a gene on both strands is kept in the annotation, and is not in $genes", {
  ann <- annotate(gtf_two_strands, "two_strands")
  expect_equal(ann$genes, ann_genc$genes[names(ann_genc$genes) != ddt])
  expect_true(ddt %in% names(ann$txs_gene))
  expect_equal(ann$trann$transcript_id, c(tx_ids, paste0(ddt_tx, "_rev")))
  # the regions of the gene have no gene
  locus <- unlist(range(ann$txs_gene[ddt]))
  regions <- ann$introns[ann$introns %over% locus]
  expect_gt(length(regions), 0)
  expect_true(all(lengths(regions$gene_id) == 0))
  # the other genes have the regions of the example
  expect_equal(ann$ncRNAs[!ann$ncRNAs %over% locus], ann_genc$ncRNAs[!ann_genc$ncRNAs %over% locus])
})

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

# Writes `lines` as a GTF file, prepares its annotation (with the other arguments of
# prepare_annotation_files() in `...`), and gives the path of the annotation file;
# `annotate()` gives the annotation
annotation_of <- function(lines, name, genome_seq = fasta, ...) {
  dir <- file.path(work, name)
  dir.create(dir)
  gtf <- file.path(dir, "genes.gtf")
  writeLines(lines, gtf)
  suppressMessages(suppressWarnings(
    prepare_annotation_files(annotation_directory = dir, gtf_file = gtf,
                             genome_seq = genome_seq, forge_BSgenome = FALSE, ...)))
  file.path(dir, "genes.gtf_Rannot")
}
annotate <- function(lines, name, ...) get(load(annotation_of(lines, name, ...)))

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

test_that("a GTF without introns gives an annotation without junctions", {
  # as a GTF of chrM, or of single-exon transcripts: each transcript gets one exon, from its transcript
  # line, and each coding transcript one CDS line, from the start of its CDS to the end of its stop codon
  type <- line_type(gtf_example)
  tx <- line_attr(gtf_example, "transcript_id")
  exons <- sub("^([^\t]*\t[^\t]*\t)transcript\t", "\\1exon\t", gtf_example[type == "transcript"])
  coding <- which(type %in% c("CDS", "stop_codon"))
  cds <- vapply(split(coding, tx[coding]), function(i) {
    x <- strsplit(gtf_example[i[type[i] == "CDS"][1]], "\t")[[1]]
    x[4:5] <- range(as.integer(unlist(lapply(strsplit(gtf_example[i], "\t"), `[`, 4:5))))
    x[8] <- "0"
    paste(x, collapse = "\t")
  }, "")
  ann <- annotate(c(gtf_example[type %in% c("gene", "transcript")], exons, cds), "no_introns")
  expect_length(ann$junctions, 0)
  expect_length(ann$introns, 0)
  expect_setequal(names(ann$cds_txs), names(cds))
})

test_that("RiboseQC_analysis() runs on a small annotation, and its P-sites profiles use the transcripts of the 5' profiles", {
  # the annotation of the GTF without introns has 140 transcripts, 49 without a CDS, and 5 selected transcripts
  # have P-sites in their CDS. Before, with fewer than 250 transcripts, the P-sites profiles used all transcripts:
  # those without a CDS, or with a short UTR, stopped the run
  annotation <- file.path(work, "no_introns", "genes.gtf_Rannot")
  skip_if_not(file.exists(annotation))
  prefix <- file.path(work, "no_introns_analysis", "chr22")
  suppressWarnings(utils::capture.output(suppressMessages(
    RiboseQC_analysis(annotation_file = annotation, bam_files = example_file("chr22_example.bam"),
                      dest_names = prefix, create_report = FALSE, rescue_all_rls = TRUE))))
  res_all <- get(load(paste0(prefix, "_results_RiboseQC_all")))
  # the profiles of the nuclear compartment exist: the run used this case
  expect_gt(length(res_all$profiles_P_sites$P_sites_bins$nucl), 0)
  expect_identical(rownames(res_all$profiles_P_sites$P_sites_bins$nucl[["all"]]),
                   rownames(res_all$profiles_fivepr$five_prime_bins$nucl[["all"]]))
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

# GENCODE's CDS lines end before the stop codon, and its stop_codon lines add it to the CDS.
# Without them, as in GTFs with CDS lines from factR, most CDS end before their stop codon.
# From ORFquant's test of this case, without its run of ORFquant
test_that("without stop_codon lines, stop_in_gtf is NA", {
  expect_identical(ann_genc$stop_in_gtf, "*")
  expect_true(is.na(annotate(gtf_example[line_type(gtf_example) != "stop_codon"], "no_stop_codons")$stop_in_gtf))
})

# The genetic code of a chromosome comes from its name: 2 for chrM, 1 for chr1 and chr22
test_that("stop_in_gtf uses the CDS on circular chromosomes if all CDS are on them", {
  expect_identical(annotate(gtf_example, "circular", circ_chroms = "chr22")$stop_in_gtf, "*")
})

test_that("stop_in_gtf uses the genetic code of the chromosome of each CDS", {
  # a genome of two copies of the region of DDT, chr1 and chrM, both not circular, and the lines of
  # DDT on each, with other ids
  lines <- gtf_example[line_attr(gtf_example, "gene_id") %in% ddt]
  fields <- strsplit(lines, "\t")
  ends <- range(as.integer(unlist(lapply(fields, `[`, 4:5))))
  ddt_on <- function(chrom) vapply(fields, function(x) {
    x[1] <- chrom
    x[4:5] <- as.integer(x[4:5]) - ends[1] + 1L
    x[9] <- gsub('(_id "[^"]*)"', paste0("\\1_", chrom, '"'), x[9])
    paste(x, collapse = "\t")
  }, "")
  region <- Rsamtools::getSeq(Rsamtools::FaFile(fasta), GenomicRanges::GRanges("chr22", IRanges::IRanges(ends[1], ends[2])))
  genome <- file.path(work, "two_codes.fa")
  Biostrings::writeXStringSet(c(stats::setNames(region, "chr1"), stats::setNames(region, "chrM")), genome)
  Rsamtools::indexFa(genome)
  ann <- annotate(c(ddt_on("chr1"), ddt_on("chrM")), "two_codes", genome_seq = genome, circ_chroms = character(0))
  expect_identical(ann$genetic_codes[c("chr1", "chrM"), "genetic_code"], c("1", "2"))
  expect_identical(ann$stop_in_gtf, "*")
  # the code of a chromosome without CDS is not used
  expect_identical(annotate(ddt_on("chr1"), "code_1", genome_seq = genome, circ_chroms = character(0))$stop_in_gtf, "*")
  # TGA, the stop codon of 5 of the 6 CDS of DDT, is W in the code of chrM
  expect_true(is.na(annotate(ddt_on("chrM"), "code_2", genome_seq = genome, circ_chroms = character(0))$stop_in_gtf))
})

test_that("load_annotation() returns the annotation invisibly, and assigns it in envir", {
  env <- new.env()
  res <- suppressMessages(withVisible(load_annotation(file.path(work, "genc", "genes.gtf_Rannot"), envir = env)))
  expect_false(res$visible)
  expect_named(res$value, c("GTF_annotation", "genome_seq"))
  expect_true(is(res$value$genome_seq, "FaFile_Circ"))
  expect_identical(res$value$GTF_annotation$trann, ann_genc$trann)
  expect_identical(env$GTF_annotation, res$value$GTF_annotation)
  expect_identical(env$genome_seq, res$value$genome_seq)
})

test_that("lines that do not agree give a warning, and the first line gives the value", {
  # a gene with two transcripts: the lines of the second get another gene name, and the last line of the
  # first another biotype
  gene_of <- line_attr(gtf_example, "gene_id")
  tx_of <- line_attr(gtf_example, "transcript_id")
  txs <- lapply(split(tx_of, gene_of), function(x) unique(stats::na.omit(x)))
  g <- names(txs)[lengths(txs) >= 2][1]
  t1 <- txs[[g]][1]
  t2 <- txs[[g]][2]
  gtf <- gtf_example
  gtf[tx_of %in% t2] <- sub('gene_name "[^"]*"', 'gene_name "OTHER_NAME"', gtf[tx_of %in% t2])
  last1 <- max(which(tx_of %in% t1))
  gtf[last1] <- sub('transcript_type "[^"]*"', 'transcript_type "OTHER_TYPE"', gtf[last1])
  dir <- file.path(work, "disagree")
  dir.create(dir)
  writeLines(gtf, file.path(dir, "genes.gtf"))
  warns <- character()
  withCallingHandlers(suppressMessages(
    prepare_annotation_files(annotation_directory = dir, gtf_file = file.path(dir, "genes.gtf"),
                             genome_seq = fasta, forge_BSgenome = FALSE)),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  warns <- grep("do not agree", warns, value = TRUE)
  expect_length(warns, 1)
  expect_match(warns, paste0("the lines of 1 transcript give different values (e.g. ", t1, " (transcript_type))"), fixed = TRUE)
  expect_match(warns, paste0("the transcripts of 1 gene give different gene names or biotypes (e.g. ", g, " (gene_name))"), fixed = TRUE)
  # the annotation is that of the first lines
  trann <- get(load(file.path(dir, "genes.gtf_Rannot")))$trann
  expect_equal(trann$transcript_biotype, ann_genc$trann$transcript_biotype)
  expect_equal(trann$gene_name[trann$transcript_id == t2], "OTHER_NAME")
  expect_equal(trann$gene_name[match(g, trann$gene_id)], ann_genc$trann$gene_name[match(g, ann_genc$trann$gene_id)])
})

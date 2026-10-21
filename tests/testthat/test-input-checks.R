# Input that cannot work must stop prepare_annotation_files() at once, with a
# message that says what is wrong. These checks need no annotation, so they run
# before any long step. From ORFquant's tests (tests/testthat/test-input-checks.R).

work <- tempfile("riboseqc_input_")
dir.create(work)
gtf_example <- readLines(gzfile(example_file("chr22_example.gtf.gz")))
fasta <- example_file("chr22_example.fa.gz")
write_gtf <- function(lines, name) {
  f <- file.path(work, name)
  writeLines(lines, f)
  f
}
# Runs prepare_annotation_files() on a GTF file and gives its error message
annotation_error <- function(gtf, genome_seq = fasta, out = file.path(work, "annotation")) {
  tryCatch({
    suppressMessages(prepare_annotation_files(annotation_directory = out, gtf_file = gtf,
                                              genome_seq = genome_seq, forge_BSgenome = FALSE))
    NA_character_
  }, error = function(e) conditionMessage(e))
}

test_that("prepare_annotation_files() names the files that do not exist, and writes nothing", {
  out <- file.path(work, "annotation_none")
  expect_match(annotation_error(file.path(work, "no.gtf"), out = out), "don't exist.*no\\.gtf")
  expect_match(annotation_error(write_gtf(gtf_example, "ok.gtf"), file.path(work, "no.fa"), out = out),
               "don't exist.*no\\.fa")
  expect_false(dir.exists(out))
})

test_that("prepare_annotation_files() asks for a genome sequence", {
  expect_error(prepare_annotation_files(annotation_directory = file.path(work, "annotation_none"),
                                        gtf_file = write_gtf(gtf_example, "ok.gtf")),
               "genome_seq")
})

test_that("prepare_annotation_files() no longer forges a BSgenome package, and writes nothing", {
  out <- file.path(work, "annotation_none")
  expect_error(prepare_annotation_files(annotation_directory = out, gtf_file = write_gtf(gtf_example, "ok.gtf"),
                                        twobit_file = file.path(work, "genome.2bit"), forge_BSgenome = TRUE),
               "genome_seq.*no longer forges a BSgenome package.*twoBitToFa")
  expect_false(dir.exists(out))
})

test_that("prepare_annotation_files() refuses a GTF file without exons, or without ids", {
  expect_match(annotation_error(write_gtf(character(), "empty.gtf")), "empty or has no exon lines")
  expect_match(annotation_error(write_gtf(gtf_example[!grepl("\texon\t", gtf_example)], "noexon.gtf")),
               "no exon lines")
  expect_match(annotation_error(write_gtf(sub('transcript_id "[^"]*"; ', "", gtf_example), "notx.gtf")),
               "no transcript_id or gene_id")
  expect_match(annotation_error(write_gtf(sub('gene_id "[^"]*"; ', "", gtf_example), "nogene.gtf")),
               "no transcript_id or gene_id")
})

test_that("prepare_annotation_files() refuses a GTF file without CDS lines, also with start and stop codons, and writes nothing", {
  # txdbmaker makes a CDS of each stop_codon line
  out <- file.path(work, "annotation_none")
  expect_match(annotation_error(write_gtf(gtf_example[!grepl("\tCDS\t", gtf_example)], "nocds.gtf"), out = out),
               "nocds\\.gtf has no CDS lines")
  expect_false(dir.exists(out))
  # the CDS lines can come before the exon lines
  cds <- grepl("\tCDS\t", gtf_example)
  expect_true(check_gtf(write_gtf(c(gtf_example[cds], gtf_example[!cds]), "cds_first.gtf")))
})

test_that("prepare_annotation_files() says when the chromosomes of the GTF are not in the genome", {
  msg <- annotation_error(write_gtf(sub("^chr22", "22", gtf_example), "nochr.gtf"))
  expect_match(msg, "Reading the GTF file .* failed")
  expect_match(msg, "chromosome names must be those of the genome sequence \\(e\\.g\\. chr22\\)")
})

test_that("check_gtf_values() warns when the lines of a transcript, or the transcripts of a gene, do not agree", {
  # a value on some lines only (as gffread writes them) agrees, and two transcripts can have two biotypes
  gtf_ids <- data.frame(gene_id = "g1", transcript_id = c(NA, "t1", "t1", "t2"),
                        transcript_type = c(NA, "protein_coding", NA, "lncRNA"))
  trann <- data.frame(gene_id = "g1", gene_biotype = "protein_coding", gene_name = "A", transcript_id = c("t1", "t2"))
  expect_silent(check_gtf_values(gtf_ids, trann))
  gtf_ids$transcript_type[3] <- "lncRNA"
  expect_warning(check_gtf_values(gtf_ids, trann),
                 "do not agree.*the lines of 1 transcript give different values \\(e\\.g\\. t1 \\(transcript_type\\)\\)")
  trann$gene_name[2] <- "B"
  msg <- tryCatch(check_gtf_values(gtf_ids, trann), warning = conditionMessage)
  expect_match(msg, "the lines of 1 transcript give")
  expect_match(msg, "the transcripts of 1 gene give different gene names or biotypes \\(e\\.g\\. g1 \\(gene_name\\)\\)")
})

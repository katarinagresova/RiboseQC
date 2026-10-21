# The genetic code of each chromosome: its name gives it, and the argument
# genetic_codes of prepare_annotation_files() replaces it.

work <- tempfile("riboseqc_codes_")
dir.create(work)
gtf_example <- readLines(gzfile(example_file("chr22_example.gtf.gz")))
fasta <- example_file("chr22_example.fa.gz")

test_that("check_genetic_codes() accepts NULL, or ids of GENETIC_CODE_TABLE named by chromosomes of the genome", {
  chroms <- c("chr1", "chrM")
  expect_silent(check_genetic_codes(NULL, chroms))
  expect_silent(check_genetic_codes(c(chrM = "3"), chroms))
  expect_silent(check_genetic_codes(c(chrM = 3, chr1 = 1), chroms))
  expect_error(check_genetic_codes("3", chroms), "one genetic code for each chromosome, named by the chromosomes")
  expect_error(check_genetic_codes(c(chrM = "3", chrM = "2"), chroms), "named by the chromosomes")
  expect_error(check_genetic_codes(list(chrM = "3"), chroms), "named by the chromosomes")
  # code 7 was deleted from the NCBI table
  expect_error(check_genetic_codes(c(chrM = "7", chr1 = "SGC0"), chroms),
               "not ids of Biostrings::GENETIC_CODE_TABLE: 7, SGC0")
  expect_error(check_genetic_codes(c(Mito = "3"), chroms),
               "not in the genome sequence: Mito \\(the genome has e\\.g\\. chr1, chrM\\)")
})

test_that("prepare_annotation_files() stops at once for genetic_codes with a chromosome that is not in the genome", {
  out <- file.path(work, "wrong")
  gtf <- file.path(work, "example.gtf")
  writeLines(gtf_example, gtf)
  expect_error(prepare_annotation_files(annotation_directory = out, gtf_file = gtf, genome_seq = fasta,
                                        genetic_codes = c(chrM = "3")),
               "genetic_codes has chromosomes that are not in the genome sequence: chrM")
  expect_false(file.exists(file.path(out, "example.gtf_TxDb")))
})

# A genome of two copies of the region of DDT, chr1 and chrM (circular), and the
# lines of DDT on each, with other ids
ddt <- "ENSG00000099977.16"
lines <- gtf_example[grepl(paste0('gene_id "', ddt, '"'), gtf_example, fixed = TRUE)]
fields <- strsplit(lines, "\t")
ends <- range(as.integer(unlist(lapply(fields, `[`, 4:5))))
ddt_on <- function(chrom) vapply(fields, function(x) {
  x[1] <- chrom
  x[4:5] <- as.integer(x[4:5]) - ends[1] + 1L
  x[9] <- gsub('(_id "[^"]*)"', paste0("\\1_", chrom, '"'), x[9])
  paste(x, collapse = "\t")
}, "")
region <- Rsamtools::getSeq(Rsamtools::FaFile(fasta), GenomicRanges::GRanges("chr22", IRanges::IRanges(ends[1], ends[2])))
genome <- file.path(work, "chr1_chrM.fa")
Biostrings::writeXStringSet(c(stats::setNames(region, "chr1"), stats::setNames(region, "chrM")), genome)
Rsamtools::indexFa(genome)
gtf <- file.path(work, "chr1_chrM.gtf")
writeLines(c(ddt_on("chr1"), ddt_on("chrM")), gtf)

# Prepares the annotation in work/name, with the other arguments in `...`, and gives it.
# txdbmaker warns that the stop_codon lines have phases, and that the TxDb has no genome
# version
annotate <- function(name, ...) {
  dir <- file.path(work, name)
  quiet(prepare_annotation_files(annotation_directory = dir, gtf_file = gtf, genome_seq = genome, ...),
        expected = c("non-NA values for features of type", "genome version information is not available"))
  get(load(file.path(dir, "chr1_chrM.gtf_Rannot")))
}

test_that("genetic_codes replaces the genetic codes that the chromosome names give, and nothing else", {
  ann <- annotate("names")
  expect_identical(ann$genetic_codes[c("chr1", "chrM"), "isCircular"], c(FALSE, TRUE))
  expect_identical(ann$genetic_codes[c("chr1", "chrM"), "genetic_code"], c("1", "2"))
  ann_yeast <- annotate("chrM_3", genetic_codes = c(chrM = 3))
  expect_identical(rownames(ann_yeast$genetic_codes), rownames(ann$genetic_codes))
  expect_identical(ann_yeast$genetic_codes[c("chr1", "chrM"), "genetic_code"], c("1", "3"))
  # the other elements are the same, apart from the creation time of the TxDb in their metadata
  no_metadata <- function(x) {
    if (methods::is(x, "Vector")) S4Vectors::metadata(x) <- list()
    x
  }
  same <- setdiff(names(ann), c("genetic_codes", "genome"))
  expect_equal(lapply(ann_yeast[same], no_metadata), lapply(ann[same], no_metadata))
  codes <- utils::read.delim(file.path(work, "chrM_3", "genetic_codes"))
  expect_identical(codes$genetic_code[codes$chromosome == "chrM"], 3L)
})

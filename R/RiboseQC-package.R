#' RiboseQC: quality control of Ribo-seq data
#'
#' RiboseQC does the quality control of Ribo-seq data and calculates the P-site offsets, for cytoplasmic
#' and organellar ribosomes (Calviello et al., bioRxiv 2019, \doi{10.1101/601468}). It works in two steps:
#' \enumerate{
#' \item \code{\link{prepare_annotation_files}} makes the annotation from a GTF file and a FASTA file of the genome.
#' \item \code{\link{RiboseQC_analysis}} analyzes the BAM files with this annotation, writes the results and makes
#' the html report (\code{\link{create_html_report}}).
#' }
#' \code{\link{prepare_for_ORFquant}} calculates the P-sites for ORFquant from a BAM file and a table of offsets.
#' The vignette runs the two steps on the example data of the package: \code{vignette("RiboseQC", package = "RiboseQC")}.
"_PACKAGE"

# Package imports. NAMESPACE is generated from these tags by roxygen2.
#' @import BSgenome
#' @rawNamespace import(BiocGenerics, except = c(combine, Position))
#' @import Biostrings
#' @import DT
#' @import GenomicAlignments
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import Rsamtools
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import knitr
#' @import reshape2
#' @import rmarkdown
#' @import rtracklayer
#' @import viridis
#' @importFrom AnnotationDbi saveDb
#' @importFrom IRanges %over% CharacterList DataFrameList IRanges IntegerList
#' @importFrom IRanges LogicalList
#' @importFrom S4Vectors runValue<- DataFrame aggregate cbind.DataFrame
#' @importFrom S4Vectors endoapply queryHits runLength runValue subjectHits
#' @importFrom Seqinfo seqlengths<- Seqinfo isCircular seqlengths seqlevels seqlevelsInUse
#' @importFrom SummarizedExperiment assay
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom methods as is new representation setMethod setRefClass
#' @importFrom stats kmeans
#' @importFrom utils read.table write.table
NULL

# Column names used in aes() and mutate(), which R CMD check takes for globals.
utils::globalVariables(c("L1", "Var1", "Var2", "aa", "biotype", "codon", "frame",
  "read_length", "value", "variable"))
# generate_rdata_list() load()s res_all.
utils::globalVariables("res_all")

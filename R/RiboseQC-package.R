# Package imports. NAMESPACE is generated from these tags by roxygen2.
#' @import BSgenome
#' @rawNamespace import(BiocGenerics, except = c(combine, Position))
#' @import Biostrings
#' @import DT
#' @import GenomicAlignments
#' @import GenomicFeatures
#' @import GenomicFiles
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
#' @importFrom Seqinfo seqlengths<- Seqinfo isCircular seqlengths seqlevels
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom methods as is new representation setMethod setRefClass
#' @importFrom stats kmeans
#' @importFrom utils write.table
NULL

# Column names used in aes() and mutate(), which R CMD check takes for globals.
utils::globalVariables(c("L1", "Var1", "Var2", "aa", "biotype", "codon", "frame",
  "read_length", "value", "variable"))
# generate_rdata_list() load()s res_all.
utils::globalVariables("res_all")

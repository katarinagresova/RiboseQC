#' @include riboseqc.R
NULL


#' Perform a Ribo-seQC analysis
#'
#' This function loads annotation created by the prepare_annotation_files function, and analyzes a BAM file.
#' @keywords Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_file Full path to the annotation file (*Rannot). Or, a vector with paths to one annotation file per bam file.
#' @param bam_files character vector containing the full path to the bam files
#' @param chunk_size the number of alignments to read at each iteration, defaults to 5000000, increase when more RAM is available. Must be more than 100000 and at most 100000000
#' @param n_cores number of cores to use, each for one chunk of \code{chunk_size} alignments at a time, in forked processes (not on Windows). Each core needs the memory for one chunk. The results do not depend on it. Defaults to 1
#' @param read_subset Calculate the offsets only for the read lengths with the most reads in each compartment: by number of reads on CDS and UTRs,
#' until they have more than 99 percent of these reads. Defaults to \code{TRUE}. Must be of length 1 or same length as bam_files.
#' @param readlength_choice_method Method used to subset relevant read lengths (see \code{choose_readlengths} function); defaults to 'max_coverage'. Must be of length 1 or same length as bam_files.
#' @param rescue_all_rls Set cutoff of 12 for all other read lengths of the bam file, in each compartment: all read lengths without a selected
#' cutoff (or without a cutoff in \code{offsets_df}). Defaults to \code{FALSE}. Must be of length 1 or same length as bam_files.
#' @param write_tmp_files Should output all the results (in *results_RiboseQC_all)? Defaults to \code{TRUE}. \code{TRUE} or \code{FALSE}, for all bam files.
#' @param dest_names character vector containing the prefixes to use for the result output files. Defaults to same as \code{bam_files}
#' @param fast_mode Use only the 500 transcripts with the most reads to build profiles (else at most 5000)? The profiles give the cutoffs and
#' the codon counts. Defaults to \code{TRUE}. Must be of length 1 or same length as bam_files.
#' @param create_report Create an html report showing the RiboseQC analysis results. Defaults to \code{TRUE}
#' @param sample_names character vector containing the names for each sample analyzed (for the html report). Defaults to 'sample1', 'sample2' ...
#' @param report_file desired filename for the html report file. Defaults to the first entry of \code{bam_files} followed by '_RiboseQC_report.html'
#' @param extended_report creates a large html report including codon occupancy for each read length. Defaults to \code{FALSE}
#' @param pdf_plots creates a pdf file for each plot of the html report, in \code{<report_file>_plots/pdf/}; only with \code{create_report = TRUE}. Defaults to \code{TRUE}
#' @param write_TSV_files write the main tables of the results as tab-separated files, for each BAM file (see \code{\link{RiboseQC_output}}). Defaults to \code{TRUE}
#' @param stranded are the analyzed libraries strand-specific? TRUE, FALSE or 'inverse'. Defaults to \code{TRUE}.
#' Use \code{'inverse'} for libraries with the reads on the strand opposite to the RNA: RiboseQC inverts the strand of each read,
#' then does the analysis of \code{TRUE}. So the P-sites and the coverage tracks are on the strand of the RNA
#' @param normalize_cov export normalized (sum to 1 million) bedgraph files for coverage tracks? Defaults to \code{TRUE}
#' @param offsets_df optionally input an offsets_df created externally. Must have cols 'read_length','cutoff'. It replaces the selected read lengths
#' and cutoffs for the P_sites, and so for \code{P_sites_stats}, \code{profiles_P_sites}, the P_sites bedgraph files, 'for_ORFquant' and 'P_sites_offsets.tsv'.
#' \code{selection_cutoffs}, \code{summary_P_sites}, 'P_sites_calcs' and the tables of read lengths and cutoffs in the report do not change: they have the calculated cutoffs.
#' A data.frame gives the offsets of the nuclear compartment (\code{"nucl"}); a list with one such table for each compartment, named as the compartments,
#' gives the offsets of the compartments with selected read lengths: it must have all of them, and the other compartments of the list are not used
#' @param genome_seq Not used: the genome sequence is the FASTA file of the annotation, the \code{genome_seq} of
#' \code{\link{prepare_annotation_files}}. A path or an \code{FaFile} of another file gives a warning. To use another FASTA
#' file, make the annotation with it
#' @return The paths of the 'results_RiboseQC' files, one for each bam file, in the order of \code{bam_files} (the \code{input_files}
#' of \code{\link{create_html_report}}). For each bam file, the function saves these files, with the prefix in \code{dest_names}
#' (by default the bam_files path):\cr
#' 'results_RiboseQC_all': an R file with the list \code{res_all} of all results. Only with \code{write_tmp_files = TRUE}.\cr
#' 'results_RiboseQC': a smaller R file used for creating the html report.\cr
#' 'P_sites_calcs': a tab-separated summary of P_sites selection statistics (\code{summary_P_sites} of \code{res_all}).\cr
#' 'P_sites_offsets.tsv': the tab-separated file with the offsets that were used to calculate the P_sites, with the columns 'read_length', 'cutoff' and 'compartment'
#' (\code{"nucl"} or the name of a circular chromosome): one row for each read length and compartment, including the read lengths of \code{offsets_df} and of \code{rescue_all_rls}.
#' It is the cutoff file that \code{\link{prepare_for_ORFquant}} reads, and it gives the P_sites_all, P_sites_uniq and P_sites_uniq_mm of 'for_ORFquant'.\cr
#' 'for_ORFquant': an R object that can be used in the ORFquant pipeline: the list \code{for_ORFquant} with the P_sites (P_sites_all, P_sites_uniq,
#' P_sites_uniq_mm) and the \code{junctions}. Only if there are P_sites.\cr
#' 'junctions': an R file with the object \code{juns}, the reads on the annotated splice junctions.\cr
#' 'P_sites_plus.bedgraph', 'P_sites_minus.bedgraph', 'P_sites_uniq_plus.bedgraph' and 'P_sites_uniq_minus.bedgraph': the P_sites positions of all and of
#' uniquely mapping reads, on each strand. Only if there are such P_sites.\cr
#' 'coverage_plus.bedgraph', 'coverage_minus.bedgraph', 'coverage_uniq_plus.bedgraph' and 'coverage_uniq_minus.bedgraph': the coverage of the reads
#' (see \code{normalize_cov}).\cr
#' With \code{create_report = TRUE}, it also saves the html report \code{report_file} and its plots (see \code{\link{create_html_report}}).
#'
#' With \code{write_TSV_files = TRUE}, also the tab-separated files 'summary.tsv', 'read_lengths.tsv', 'regions.tsv', 'genes.tsv',
#' 'P_sites_calcs.tsv' and 'top50_positions.tsv', with the main tables of \code{res_all}.
#'
#' \code{\link{RiboseQC_output}} describes each file and each element of \code{res_all}.
#' @details This function loads different genomic regions created in the \code{prepare_annotation_files} step,
#' separating features on different recognized organelles. The bam files is then analyzed in chunks to minimize RAM usage.\cr
#' The results of each bam file are in the list \code{res_all}, described in \code{\link{RiboseQC_output}}.
#'
#' Alignments flagged as duplicates, secondary alignments and alignments with an insertion or a deletion (\code{I} or \code{D} in the CIGAR string)
#' are not used: the read counts, the read length distributions, the offsets, the P_sites, the coverage and the junction reads do not include them,
#' and the report does not count them.
#'
#' Soft-clipped bases (\code{S} in the CIGAR string) are not counted: the read length is the length of the aligned part of the read,
#' and the offset (cutoff) counts from its first aligned base at the 5' end. So a 31-nt read with 2 soft-clipped bases counts as a 29-nt read.
#' Aligners that soft-clip reads, for example STAR with its default \code{--alignEndsType Local}, can give other read lengths, offsets and P_sites
#' than an end-to-end alignment of the same reads (STAR: \code{--alignEndsType EndToEnd}). Remove adapters and UMIs before the alignment.
#' An \code{offsets_df} must come from the same alignments.
#' @seealso \code{\link{RiboseQC_output}}, \code{\link{prepare_annotation_files}}, \code{\link{calc_cutoffs_from_profiles}}, \code{\link{choose_readlengths}}, \code{\link{create_html_report}}.
#' @export 

RiboseQC_analysis <- function(annotation_file, bam_files, read_subset = TRUE, readlength_choice_method = "max_coverage", 
    genome_seq = NULL, stranded = TRUE, normalize_cov = TRUE, chunk_size = 5000000L, 
    write_tmp_files = TRUE, dest_names = NA, rescue_all_rls = FALSE, fast_mode = TRUE, 
    create_report = TRUE, sample_names = NA, report_file = NA, extended_report = FALSE, 
    pdf_plots = TRUE, offsets_df = NULL, n_cores = 1, write_TSV_files = TRUE) {

    check_n_cores(n_cores)

    if (length(dest_names) == 1) {
        if (is.na(dest_names)) {
            dest_names = bam_files
        }
    }
    
    if (sum(duplicated(dest_names)) > 0) {
        stop(paste("duplicated 'dest_names' parameter (possibly same bam file). Please specify different 'dest_names'.", 
            date(), "\n"))
    }
    
    if (is.na(report_file)) {
        report_file = paste(bam_files[1], "_RiboseQC_report.html", sep = "")
    }
    if (length(sample_names) == 1) {
        if (is.na(sample_names)) {
            sample_names = paste("sample", seq_along(bam_files), sep = "")
        }
    }
    
    if (!is.logical(read_subset)) {
        stop(paste("'read_subset' must be logical (TRUE or FALSE, no quotes).", date(), 
            "\n"))
    }
    if (!is.logical(write_tmp_files)) {
        stop(paste("'write_tmp_files' must be logical (TRUE or FALSE, no quotes).", 
            date(), "\n"))
    }
    if (!is.logical(rescue_all_rls)) {
        stop(paste("'rescue_all_rls' must be logical (TRUE or FALSE, no quotes).", 
            date(), "\n"))
    }
    if (!is.logical(fast_mode)) {
        stop(paste("'fast_mode' must be logical (TRUE or FALSE, no quotes).", date(), 
            "\n"))
    }
    if (!is.logical(create_report)) {
        stop(paste("'create_report' must be logical (TRUE or FALSE, no quotes).", 
            date(), "\n"))
    }
    # one value: else the run stops when it writes the files of the first BAM file
    if (!isTRUE(write_TSV_files) && !isFALSE(write_TSV_files)) {
        stop(paste("'write_TSV_files' must be TRUE or FALSE (no quotes).", date(), "\n"))
    }
    if (!is.integer(chunk_size) & !is.numeric(chunk_size)) {
        stop(paste("'chunk_size' must be an integer value (e.g. 5000000).", date(), 
            "\n"))
    }
    if (chunk_size <= 1e+05 | chunk_size > 1e+08) {
        stop(paste("'chunk_size' must be an integer value between 100000 and 100000000 (e.g. 5000000).", 
            date(), "\n"))
    }
    if (!is.character(annotation_file)) {
        stop(paste("'annotation_file' must be a character (e.g.'/home/myfolder/myfile').", 
            date(), "\n"))
    }
    if (!is.character(bam_files)) {
        stop(paste("'bam_files' must be a character (e.g.'/home/myfolder/myfile') or character vector (e.g. c('/home/myfolder/myfile1','/home/myfolder/myfile2').", 
            date(), "\n"))
    }
    if (!is.character(dest_names)) {
        stop(paste("'dest_names' must be a character (e.g.'/home/myfolder/myfile') or character vector (e.g. c('/home/myfolder/myfile1','/home/myfolder/myfile2')..", 
            date(), "\n"))
    }
    if (!is.character(report_file)) {
        stop(paste("'report_file' must be a character (e.g.'/home/myfolder/myfile').", 
            date(), "\n"))
    }
    if (!is.character(sample_names)) {
        stop(paste("'sample_names' must be a character (e.g.'failed_experiment'). or character vector (e.g. c('attempt_1','attempt_2').", 
            date(), "\n"))
    }
    if (!is.logical(normalize_cov)) {
        stop(paste("'normalize_cov' must be logical (TRUE or FALSE, no quotes).", 
            date(), "\n"))
    }
    
    if (!length(read_subset) %in% c(1, length(bam_files))) {
        stop(paste("length of 'read_subset' must be 1 or same as 'bam_files'.", date(), 
            "\n"))
    }
    if (!length(readlength_choice_method) %in% c(1, length(bam_files))) {
        stop(paste("length of 'readlength_choice_method' must be 1 or same as 'bam_files'.", 
            date(), "\n"))
    }
    if (!length(rescue_all_rls) %in% c(1, length(bam_files))) {
        stop(paste("length of 'rescue_all_rls' must be 1 or same as 'bam_files'.", 
            date(), "\n"))
    }
    if (!length(fast_mode) %in% c(1, length(bam_files))) {
        stop(paste("length of 'fast_mode' must be 1 or same as 'bam_files'.", date(), 
            "\n"))
    }
    if (!length(stranded) %in% c(1, length(bam_files))) {
        stop(paste("length of 'stranded' must be 1 or same as 'bam_files'.", date(), 
            "\n"))
    }
    if (!length(normalize_cov) %in% c(1, length(bam_files))) {
        stop(paste("length of 'normalize_cov' must be 1 or same as 'bam_files'.", 
            date(), "\n"))
    }
    if (length(bam_files) > 1) {
        if (length(read_subset) == 1) {
            read_subset <- rep(read_subset, length(bam_files))
        }
        if (length(readlength_choice_method) == 1) {
            readlength_choice_method <- rep(readlength_choice_method, length(bam_files))
        }
        if (length(rescue_all_rls) == 1) {
            rescue_all_rls <- rep(rescue_all_rls, length(bam_files))
        }
        if (length(fast_mode) == 1) {
            fast_mode <- rep(fast_mode, length(bam_files))
        }
        if (length(stranded) == 1) {
            stranded <- rep(stranded, length(bam_files))
        }
        if (length(normalize_cov) == 1) {
            normalize_cov <- rep(normalize_cov, length(bam_files))
        }
    }
    
    
    
    if (!length(dest_names) %in% c(length(bam_files))) {
        stop(paste("length of 'dest_names' must be same as 'bam_files'.", date(), 
            "\n"))
    }
    if (!length(sample_names) %in% c(length(bam_files))) {
        stop(paste("length of 'sample_names' must be same as 'bam_files'.", date(), 
            "\n"))
    }
    if (!length(annotation_file) %in% c(1, length(bam_files))) {
        stop(paste("length of 'annotation_file' must be 1 or same as 'bam_files'.", 
            date(), "\n"))
    }
    
    if (sum(!readlength_choice_method %in% c("max_coverage", "max_inframe", "all")) > 
        0) {
        stop(paste("'readlength_choice_method' must be one of 'max_coverage','max_inframe','all'.", 
            date(), "\n"))
    }
    # stranded is "inverse", or a value that as.logical() reads (TRUE, FALSE, "TRUE", ...): else the run stopped
    # in the counting pass of the first BAM file, after the pass of the read lengths over all its reads
    bad_stranded <- !stranded %in% "inverse" & is.na(as.logical(stranded))
    if (any(bad_stranded)) {
        stop(paste("'stranded' must be TRUE, FALSE or 'inverse', not ", paste(unique(stranded[bad_stranded]),
            collapse = ", "), ". ", date(), "\n", sep = ""))
    }
    
    # all input files, before the first BAM file: else a missing file stops the run after the
    # analysis of the BAM files before it
    check_files_exist(c(annotation_file, bam_files))
    
    # genome_seq is not used: the genome sequence is the FASTA file of the annotation, which load_annotation()
    # assigns to genome_seq. A warning says so when genome_seq is another file
    given_genome_seq <- genome_seq
    
    list_annotations <- list()
    
    for (annots in seq_along(annotation_file)) {
        message(paste("Loading annotation files in ", annotation_file[annots], " ... ", 
            date(), sep = ""))
        
        load_annotation(annotation_file[annots])
        check_genome_seq(given_genome_seq, genome_seq, annotation_file[annots])
        # the annotation refers to its FASTA file by its path, and the P-sites profiles, near the end of
        # the analysis of each BAM file, are the first step that reads it
        if (is(genome_seq, "FaFile") && !file.exists(genome_seq$path)) {
            stop(paste("The FASTA file of the annotation ", annotation_file[annots], " doesn't exist: ", 
                genome_seq$path, ". The annotation refers to it by this path: move the file back, or make the annotation again. ", 
                date(), "\n", sep = ""))
        }
        lst <- list(GTF_annotation, genome_seq)
        names(lst) <- c("GTF_annotation", "genome_seq")
        list_annotations[[annots]] <- lst
        
    }
    
    message(paste("Loading annotation files --- Done! ", date(), sep = ""))
    
    # the chromosomes of each BAM file and of its annotation, before the first BAM file
    for (bammo in seq_along(bam_files)) {
        annots <- if (length(annotation_file) > 1) bammo else 1
        check_bam_chromosomes(seqinfo(BamFile(bam_files[bammo])), list_annotations[[annots]]$GTF_annotation$seqinfo, 
            bam_files[bammo], annotation_file[annots])
        warn_bam_chromosomes(bam_files[bammo], list_annotations[[annots]]$GTF_annotation$seqinfo, annotation_file[annots])
    }
    
    # the output directories, which are made if they don't exist, as before, must be writable: else
    # the run stops at the first file that it writes, and at the report after the analysis
    out_files <- c(dest_names, if (create_report) report_file)
    for (out_dir in unique(dirname(out_files))) {
        suppressWarnings(dir.create(out_dir, recursive = TRUE))
    }
    for (out_file in out_files) {
        check_output_dir(out_file)
    }
    
    GTF_annotation <- list_annotations[[1]]$GTF_annotation
    genome_seq <- list_annotations[[1]]$genome_seq
    
    
    resfilelist <- c()
    
    for (bammo in seq_along(bam_files)) {
        
        chunk_size <- as.integer(chunk_size)
        
        GTF_annotation <- list_annotations[[1]]$GTF_annotation
        genome_seq <- list_annotations[[1]]$genome_seq
        
        if (length(annotation_file) > 1 & bammo > 1) {
            GTF_annotation <- list_annotations[[bammo]]$GTF_annotation
            genome_seq <- list_annotations[[bammo]]$genome_seq
        }
        
        bam_file <- bam_files[bammo]
        dest_name <- dest_names[bammo]
        
        opts <- BamFile(file = bam_file, yieldSize = chunk_size)  #read BAMfile in one chunk
        # without the MD tag: only the P-sites use it, with their own param
        param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE), 
            what = c("mapq"))
        
        dira <- paste(dirname(dest_name), paste("tmp_RiboseQC", basename(dest_name), 
            sep = "_"), sep = "/")
        suppressWarnings(dir.create(dira, recursive = TRUE))
        
        seqs <- seqinfo(opts)
        circs_seq <- seqnames(GTF_annotation$seqinfo)[which(isCircular(GTF_annotation$seqinfo))]
        circs <- seqs@seqnames[which(seqs@seqnames %in% circs_seq)]
        
        red_cdss <- reduce(GTF_annotation$cds_genes)
        
        red_ex <- unlist(GTF_annotation$exons_txs)
        red_ex$gene_id <- GTF_annotation$trann$gene_id[match(names(red_ex), GTF_annotation$trann$transcript_id)]
        red_ex <- reduce(split(red_ex, red_ex$gene_id))
        
        regions <- list(reduce(unlist(GTF_annotation$cds_genes)), GTF_annotation$fiveutrs, 
            GTF_annotation$threeutrs, GTF_annotation$ncIsof, GTF_annotation$ncRNAs, 
            GTF_annotation$introns, GTF_annotation$intergenicRegions)
        names(regions) <- c("cds", "fiveutrs", "threeutrs", "ncIsof", "ncRNAs", "introns", 
            "intergenic")
        
        regions <- GRangesList(endoapply(regions, function(x) {
            mcols(x) <- NULL
            x
        }))
        
        list_locat <- list()
        list_locat[["nucl"]] <- regions
        all <- regions
        if (length(circs) > 0) {
            for (i in c("nucl", circs)) {
                if (i == "nucl") {
                  regions <- all
                  for (w in seq_along(all)) {
                    regions[[w]] <- regions[[w]][!seqnames(regions[[w]]) %in% circs]
                  }
                  list_locat[["nucl"]] <- regions
                }
                if (i != "nucl") {
                  regions <- all
                  for (w in seq_along(all)) {
                    regions[[w]] <- regions[[w]][seqnames(regions[[w]]) %in% i]
                  }
                  list_locat[[i]] <- regions
                }
            }
        }
        
        # get readlengths
        
        reduc <- function(x, y) {
            list(max(x[[1]], y[[1]]), min(x[[2]], y[[2]]))
        }
        # the read lengths need only the cigars: the flags of param, without its other fields
        yiel <- function(x) {
            readGAlignments(x, param = ScanBamParam(flag = bamFlag(param, asInteger = TRUE)))
        }
        
        mapp <- function(x) {
            # the reads with I or D go: one grepl(), not two grep() for each
            indel <- grepl("[ID]", cigar(x))
            if (any(indel)) {
                x <- x[!indel]
            }
            
            
            # the read length without the soft clips
            len_adj <- cigarillo::cigar_extent_along_query(x@cigar, after.soft.clipping = TRUE)
            mcols(x)$len_adj <- len_adj
            
            
            maxr <- max(mcols(x)$len_adj)
            minr <- min(mcols(x)$len_adj)
            list(maxr, minr)
        }
        
        message(paste("Extracting read lengths from ", bam_file, " ... ", date(), 
            sep = ""))
        maxmin <- reduce_by_yield(X = opts, YIELD = yiel, MAP = mapp, REDUCE = reduc, n_cores = n_cores)
        message(paste("Extracting read lengths --- Done! ", date(), sep = ""))
        # no chunk: else the next line stops with "subscript out of bounds"
        if (length(maxmin) == 0) {
            stop(paste("The BAM file ", bam_file, " has no mapped reads (secondary alignments and duplicates are not used). ", 
                date(), "\n", sep = ""))
        }
        readlengths <- seq(maxmin[[2]], maxmin[[1]])
        
        
        
        # body of the iteration for reducebyYield
        strandedness <- stranded[bammo]
        
        # what to do with chunks: x present chunk, y old chunks (cumulative)
        reduc <- function(x, y) {
            # the 5' ends of the chunks are summed after the last chunk, all at once (sum_chunks()): the
            # merge of each chunk with the sums of the chunks before it took longer with each chunk
            chunks_pos1 <- x[["chunks_pos1"]]
            if (is.null(chunks_pos1)) {
                chunks_pos1 <- list(x[["reads_pos1"]])
            }
            chunks_pos1 <- c(chunks_pos1, list(y[["reads_pos1"]]))
            
            
            cntsss <- y[["counts_cds_genes"]]
            cntsss$reads <- cntsss$reads + x[["counts_cds_genes"]]$reads
            
            cntsss_all <- y[["counts_all_genes"]]
            cntsss_all$reads <- cntsss_all$reads + x[["counts_all_genes"]]$reads
            
            
            cntsss_unq <- y[["counts_cds_genes_unq"]]
            cntsss_unq$reads <- cntsss_unq$reads + x[["counts_cds_genes_unq"]]$reads
            
            cntsss_all_unq <- y[["counts_all_genes_unq"]]
            cntsss_all_unq$reads <- cntsss_all_unq$reads + x[["counts_all_genes_unq"]]$reads
            
            
            reads_summary <- mapply(FUN = function(x, y) {
                DataFrame(as.matrix(x) + as.matrix(y))
            }, x[["reads_summary"]], y[["reads_summary"]], SIMPLIFY = FALSE)
            reads_summary_unq <- mapply(FUN = function(x, y) {
                DataFrame(as.matrix(x) + as.matrix(y))
            }, x[["reads_summary_unq"]], y[["reads_summary_unq"]], SIMPLIFY = FALSE)
            
            lis <- list((x[["rld"]] + y[["rld"]]), x[["rld_unq"]] + y[["rld_unq"]], 
                (x[["positions"]] + y[["positions"]]), (x[["positions_unq"]] + y[["positions_unq"]]), 
                x[["reads_pos1"]], cntsss, cntsss_unq, cntsss_all, cntsss_all_unq, reads_summary, 
                reads_summary_unq, chunks_pos1)
            names(lis) <- c("rld", "rld_unq", "positions", "positions_unq", "reads_pos1", 
                "counts_cds_genes", "counts_cds_genes_unq", "counts_all_genes", "counts_all_genes_unq", 
                "reads_summary", "reads_summary_unq", "chunks_pos1")
            lis
        }
        
        # what to do with each chunk (read as alignment file)
        yiel <- function(x) {
            readGAlignments(x, param = param)
        }
        
        # operations on the chunk (here count reads and whatnot)
        mapp <- function(x) {
            if (strandedness == "inverse") {
                x <- invertStrand(x)
                # the inverted reads are on the strand of the RNA: count them as a stranded library
                strandedness <- T
            }
            
            # Remove Insertion and Deletion (TBA)
            
            # the reads with I or D go: one grepl(), not two grep() for each
            indel <- grepl("[ID]", cigar(x))
            if (any(indel)) {
                x <- x[!indel]
            }
            
            emptyt <- rep(0, length(readlengths))
            names(emptyt) <- paste("reads", readlengths, sep = "_")
            
            
            # softclipping part
            
            
            # the read length without the soft clips
            len_adj <- cigarillo::cigar_extent_along_query(x@cigar, after.soft.clipping = TRUE)
            mcols(x)$len_adj <- len_adj
            
            # Remove S from Cigar (read positions/length are already adjusted) it helps
            # calculating P-sites positions for spliced reads
            
            cigg <- cigar(x)
            cigg_s <- grep(cigg, pattern = "S")
            if (length(cigg_s) > 0) {
                cigs <- cigg[cigg_s]
                cigs <- gsub(cigs, pattern = "^[0-9]+S", replacement = "")
                cigs <- gsub(cigs, pattern = "[0-9]+S$", replacement = "")
                cigg[cigg_s] <- cigs
                x@cigar <- cigg
            }
            mcols(x)$cigar_str <- x@cigar
            x_uniq <- x[x@elementMetadata$mapq > 50]
            
            
            
            # rld_len, rld_loc
            
            
            # initialise read length vector per compartment
            list_vects <- list(nucl = emptyt)
            if (length(circs) > 0) {
                for (i in circs) {
                  list_vects[[i]] <- emptyt
                }
            }
            
            list_vects_unq <- list(nucl = emptyt)
            if (length(circs) > 0) {
                for (i in circs) {
                  list_vects_unq[[i]] <- emptyt
                }
            }
            
            # sort reads per compartment
            list_reads <- list(nucl = x)
            list_reads_unq <- list(nucl = x_uniq)
            if (length(circs) > 0) {
                nucl <- subset(x, !seqnames(x) %in% circs)
                nucl_unq <- subset(x_uniq, !seqnames(x_uniq) %in% circs)
                list_reads <- list(nucl = nucl)
                list_reads_unq <- list(nucl = nucl_unq)
                for (i in circs) {
                  list_reads[[i]] <- subset(x, seqnames(x) == i)
                  list_reads_unq[[i]] <- subset(x_uniq, seqnames(x_uniq) == i)
                }
            }
            
            # count reads per read length per compartment
            
            for (i in names(list_reads)) {
                tab <- table(mcols(list_reads[[i]])$len_adj)
                tab_unq <- table(mcols(list_reads_unq[[i]])$len_adj)
                for (j in readlengths) {
                  cont <- as.numeric(tab[which(names(tab) == j)])
                  cont_unq <- as.numeric(tab_unq[which(names(tab_unq) == j)])
                  if (length(cont) > 0) {
                    list_vects[[i]][j - (min(readlengths) - 1)] <- cont
                  }
                  if (length(cont_unq) > 0) {
                    list_vects_unq[[i]][j - (min(readlengths) - 1)] <- cont_unq
                  }
                  
                }
            }
            
            rld_len <- do.call(what = rbind.data.frame, list_vects)
            names(rld_len) <- paste("reads", readlengths, sep = "_")
            row.names(rld_len) <- names(list_reads)
            
            rld_len_unq <- do.call(what = rbind.data.frame, list_vects_unq)
            names(rld_len_unq) <- paste("reads", readlengths, sep = "_")
            row.names(rld_len_unq) <- names(list_reads)
            
            
            # the counts of summarizeOverlaps(mode = "Union", inter.feature = FALSE) for the reads x and for the
            # unique reads x_uniq, from one findOverlaps() with x: x_uniq are the reads of x with mapq > 50
            count_reads <- function(features) {
                hits <- findOverlaps(features, x, ignore.strand = !as.logical(strandedness))
                unq <- x@elementMetadata$mapq[subjectHits(hits)] > 50
                list(all = matrix(tabulate(queryHits(hits), length(features)), dimnames = list(names(features), "reads")),
                  unq = matrix(tabulate(queryHits(hits)[unq], length(features)), dimnames = list(names(features), "reads")))
            }
            
            # count reads per biotype location per compartment
            list_rld_loc <- list()
            list_rld_loc_unq <- list()
            for (i in c("nucl", circs)) {
                cnts <- count_reads(GRangesList(list_locat[[i]]))
                list_rld_loc[[i]] <- cnts$all
                list_rld_loc_unq[[i]] <- cnts$unq
            }
            rld_loc <- do.call(what = cbind.data.frame, list_rld_loc)
            names(rld_loc) <- paste("reads", names(list_rld_loc), sep = "_")
            
            
            rld_loc_unq <- do.call(what = cbind.data.frame, list_rld_loc_unq)
            names(rld_loc_unq) <- paste("reads", names(list_rld_loc_unq), sep = "_")
            
            
            # reads_summary
            
            reads_summary <- list()
            
            # iterate over genome types and compartments
            for (i in names(list_reads)) {
                reads_by_length <- list()
                rdss <- list_reads[[i]]
                
                rgs <- GRangesList(list_locat[[i]])
                # one findOverlaps() for the reads of all read lengths, then the counts of each read length: the
                # counts of summarizeOverlaps(mode = "Union", inter.feature = FALSE) for each read length
                ovs <- suppressWarnings(findOverlaps(rgs, rdss, ignore.strand = FALSE))
                ovs <- lapply(split(queryHits(ovs), factor(mcols(rdss)$len_adj)[subjectHits(ovs)]), FUN = function(x) {
                  matrix(tabulate(x, length(rgs)), dimnames = list(names(rgs), "reads"))
                })
                nope <- readlengths[which(!readlengths %in% names(ovs))]
                if (length(nope) > 0) {
                  for (j in nope) {
                    nop <- matrix(0, nrow = 7, ncol = 1)
                    rownames(nop) <- names(rgs)
                    ovs[[as.character(j)]] <- nop
                    
                  }
                }
                ovs <- ovs[names(ovs) %in% as.character(readlengths)]
                ovs <- ovs[as.character(readlengths)]
                ovs <- do.call(ovs, what = cbind)
                colnames(ovs) <- paste("reads", readlengths, sep = "_")
                reads_summary[[i]] <- DataFrame(ovs)
                
                
            }
            
            
            
            reads_summary_unq <- list()
            
            # iterate over genome types and compartments
            for (i in names(list_reads_unq)) {
                reads_by_length <- list()
                rdss <- list_reads_unq[[i]]
                
                rgs <- GRangesList(list_locat[[i]])
                # one findOverlaps() for the reads of all read lengths, then the counts of each read length: the
                # counts of summarizeOverlaps(mode = "Union", inter.feature = FALSE) for each read length
                ovs <- suppressWarnings(findOverlaps(rgs, rdss, ignore.strand = FALSE))
                ovs <- lapply(split(queryHits(ovs), factor(mcols(rdss)$len_adj)[subjectHits(ovs)]), FUN = function(x) {
                  matrix(tabulate(x, length(rgs)), dimnames = list(names(rgs), "reads"))
                })
                nope <- readlengths[which(!readlengths %in% names(ovs))]
                if (length(nope) > 0) {
                  for (j in nope) {
                    nop <- matrix(0, nrow = 7, ncol = 1)
                    rownames(nop) <- names(rgs)
                    ovs[[as.character(j)]] <- nop
                    
                  }
                }
                ovs <- ovs[names(ovs) %in% as.character(readlengths)]
                ovs <- ovs[as.character(readlengths)]
                ovs <- do.call(ovs, what = cbind)
                colnames(ovs) <- paste("reads", readlengths, sep = "_")
                reads_summary_unq[[i]] <- DataFrame(ovs)
                
                
            }
            
            
            # reads_pos1_stst
            
            
            # take first position of each read (without the read length as the name of each range, as after the
            # merge of 2 or more chunks)
            reads_pos1 <- unlist(resize(split(GRanges(x), f = mcols(x)$len_adj), 
                1), use.names = FALSE)
            
            reads_pos1 <- split(reads_pos1, reads_pos1$len_adj)
            
            reads_pos1 <- GRangesList(lapply(reads_pos1, function(y) {
                unq <- unique(y)
                mcols(unq) <- NULL
                unq$score <- countOverlaps(unq, y, type = "equal")
                unq
                
            }))
            
            # cnts_cds_genes, cnts_all_genes
            
            seqgene <- cbind(as.vector(seqnames(GTF_annotation$genes)), GTF_annotation$genes$gene_id)
            
            
            cnts_cds <- count_reads(red_cdss)
            cnts_all <- count_reads(red_ex)
            
            cnts_cds_genes <- cnts_cds$all
            
            cnts_all_genes <- cnts_all$all
            
            cnts_cds_genes_unq <- cnts_cds$unq
            
            cnts_all_genes_unq <- cnts_all$unq
            
            
            chrsss <- seqgene[match(rownames(cnts_cds_genes), seqgene[, 2]), 1]
            
            cnts_cds_genes <- DataFrame(cnts_cds_genes)
            cnts_cds_genes$chromosome <- chrsss
            cnts_cds_genes$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_cds_genes), 
                GTF_annotation$trann$gene_id)]
            cnts_cds_genes$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_cds_genes), 
                GTF_annotation$trann$gene_id)]
            
            cnts_all_genes <- DataFrame(cnts_all_genes)
            cnts_all_genes$chromosome <- seqgene[match(rownames(cnts_all_genes), 
                seqgene[, 2]), 1]
            cnts_all_genes$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_all_genes), 
                GTF_annotation$trann$gene_id)]
            cnts_all_genes$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_all_genes), 
                GTF_annotation$trann$gene_id)]
            
            
            cnts_cds_genes_unq <- DataFrame(cnts_cds_genes_unq)
            cnts_cds_genes_unq$chromosome <- chrsss
            cnts_cds_genes_unq$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_cds_genes_unq), 
                GTF_annotation$trann$gene_id)]
            cnts_cds_genes_unq$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_cds_genes_unq), 
                GTF_annotation$trann$gene_id)]
            
            cnts_all_genes_unq <- DataFrame(cnts_all_genes_unq)
            cnts_all_genes_unq$chromosome <- seqgene[match(rownames(cnts_all_genes_unq), 
                seqgene[, 2]), 1]
            cnts_all_genes_unq$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_all_genes_unq), 
                GTF_annotation$trann$gene_id)]
            cnts_all_genes_unq$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_all_genes_unq), 
                GTF_annotation$trann$gene_id)]
            
            
            
            chunk_res <- list(rld_len, rld_len_unq, rld_loc, rld_loc_unq, reads_pos1, 
                cnts_cds_genes, cnts_cds_genes_unq, cnts_all_genes, cnts_all_genes_unq, 
                reads_summary, reads_summary_unq)
            # export list to use in the reduc part above
            names(chunk_res) <- c("rld", "rld_unq", "positions", "positions_unq", 
                "reads_pos1", "counts_cds_genes", "counts_cds_genes_unq", "counts_all_genes", 
                "counts_all_genes_unq", "reads_summary", "reads_summary_unq")
            
            chunk_res
        }
        
        
        # main for reducebyYield
        
        message(paste("Analyzing BAM file:", bam_file, "...", date()))
        
        read_stats <- reduce_by_yield(X = opts, YIELD = yiel, MAP = mapp, REDUCE = reduc,
            n_cores = n_cores)
        # two or more chunks: sum the 5' ends of all chunks
        if (!is.null(read_stats[["chunks_pos1"]])) {
            seql <- seqlevels(GTF_annotation$seqinfo)
            seqle <- seqlengths(GTF_annotation$seqinfo)
            # and the chromosomes with reads that the annotation does not have, with their lengths in the
            # BAM file: one chunk keeps these reads, so more chunks must keep them too
            seq_extra <- unlist(lapply(read_stats$chunks_pos1, seqlevelsInUse))
            seq_extra <- seqlevels(seqs)[seqlevels(seqs) %in% setdiff(seq_extra, seql)]
            seql <- c(seql, seq_extra)
            seqle <- c(seqle, seqlengths(seqs)[seq_extra])
            read_stats$reads_pos1 <- sum_chunks(read_stats$chunks_pos1, seql, seqle)
            read_stats$chunks_pos1 <- NULL
        }
        names(read_stats) <- c("rld", "rld_unq", "positions", "positions_unq", "reads_pos1", 
            "counts_cds_genes", "counts_cds_genes_unq", "counts_all_genes", "counts_all_genes_unq", 
            "reads_summary", "reads_summary_unq")
        
        cds_cnts <- read_stats$counts_cds_genes
        cds_cnts$RPKM <- cds_cnts$reads/(sum(cds_cnts$reads)/1e+06)
        cds_cnts$RPKM <- round(cds_cnts$RPKM/(sum(width(red_cdss))/1000), digits = 4)
        cds_cnts$TPM <- cds_cnts$reads/(sum(width(red_cdss))/1000)
        cds_cnts$TPM <- round(cds_cnts$TPM/(sum(cds_cnts$TPM)/1e+06), digits = 4)
        read_stats$counts_cds_genes <- cds_cnts
        
        cds_cnts_unq <- read_stats$counts_cds_genes_unq
        cds_cnts_unq$RPKM <- cds_cnts_unq$reads/(sum(cds_cnts_unq$reads)/1e+06)
        cds_cnts_unq$RPKM <- round(cds_cnts_unq$RPKM/(sum(width(red_cdss))/1000), 
            digits = 4)
        cds_cnts_unq$TPM <- cds_cnts_unq$reads/(sum(width(red_cdss))/1000)
        cds_cnts_unq$TPM <- round(cds_cnts_unq$TPM/(sum(cds_cnts_unq$TPM)/1e+06), 
            digits = 4)
        read_stats$counts_cds_genes_unq <- cds_cnts_unq
        
        ex_cnts <- read_stats$counts_all_genes
        ex_cnts$RPKM <- ex_cnts$reads/(sum(ex_cnts$reads)/1e+06)
        ex_cnts$RPKM <- round(ex_cnts$RPKM/(sum(width(red_ex))/1000), digits = 4)
        ex_cnts$TPM <- ex_cnts$reads/(sum(width(red_ex))/1000)
        ex_cnts$TPM <- round(ex_cnts$TPM/(sum(ex_cnts$TPM)/1e+06), digits = 4)
        read_stats$counts_all_genes <- ex_cnts
        
        ex_cnts_unq <- read_stats$counts_all_genes_unq
        ex_cnts_unq$RPKM <- ex_cnts_unq$reads/(sum(ex_cnts_unq$reads)/1e+06)
        ex_cnts_unq$RPKM <- round(ex_cnts_unq$RPKM/(sum(width(red_ex))/1000), digits = 4)
        ex_cnts_unq$TPM <- ex_cnts_unq$reads/(sum(width(red_ex))/1000)
        ex_cnts_unq$TPM <- round(ex_cnts_unq$TPM/(sum(ex_cnts_unq$TPM)/1e+06), digits = 4)
        read_stats$counts_all_genes_unq <- ex_cnts_unq
        
        save(read_stats, file = paste(dira, "read_stats", sep = "/"))
        message(paste("Analyzing BAM file --- Done!", date()))
        
        ps <- read_stats$reads_pos1
        
        tile_cds <- GTF_annotation$cds_txs_coords
        tile_cds <- tile_cds[tile_cds$reprentative_mostcommon]
        
        list_txs_ok <- GRangesList()
        for (ci in circs) {
            citxs <- GTF_annotation$cds_txs[seqnames(GTF_annotation$cds_txs) == ci]
            citxs <- citxs[elementNROWS(citxs) > 0]
            citxs_txs <- GTF_annotation$cds_txs_coords[as.vector(seqnames(GTF_annotation$cds_txs_coords)) %in% 
                names(citxs)]
            list_txs_ok[[ci]] <- citxs_txs
        }
        list_txs_ok[["nucl"]] <- tile_cds[!tile_cds %in% unlist(list_txs_ok)]
        
        message(paste("Building aggregate 5' profiles ...", date()))
        ps_signals_win_all <- list()
        ps_signals_tiles_all <- list()
        
        # calculate profiles for bins and windows (single nt resolution)
        for (comp in c("nucl", circs)) {
            tile_cds <- list_txs_ok[[comp]]
            strand(tile_cds) <- "+"
            tile_cds <- tile_cds[width(tile_cds) > 100]
            checkgen <- TRUE
            if (comp == "nucl") {
                if ((sum((start(tile_cds) > 50))/length(start(tile_cds))) > 0.33) {
                  tile_cds <- tile_cds[start(tile_cds) > 50]
                  checkgen <- FALSE
                }
                if (sum(seqlengths(tile_cds)[as.vector(seqnames(tile_cds))] - end(tile_cds) > 
                  100)/length(tile_cds) > 0.33) {
                  tile_cds <- tile_cds[seqlengths(tile_cds)[as.vector(seqnames(tile_cds))] - 
                    end(tile_cds) > 100]
                  checkgen <- FALSE
                }
                ps_comp <- GRangesList(lapply(ps, function(x) {
                  sort(x[!seqnames(x) %in% circs])
                }))
            } else {
                ps_comp <- GRangesList(lapply(ps, function(x) {
                  sort(x[seqnames(x) == comp])
                }))
            }
            
            # select only some read lenghts, otherwise too much output and computation
            
            summm <- read_stats$reads_summary
            a <- unlist(summm[[comp]]["cds", ]) + unlist(summm[[comp]]["fiveutrs", 
                ]) + unlist(summm[[comp]]["threeutrs", ])
            names(a) <- gsub(names(a), pattern = "reads_", replacement = "")
            a <- sort(a/(sum(a)/100), decreasing = TRUE)
            cuma <- cumsum(a)
            if (length(cuma) == 0) {
                ps_signals_tiles_all[[comp]] <- c()
                ps_signals_win_all[[comp]] <- c()
                next
            }
            # if read subset, discard low abundance
            if (read_subset[bammo] == TRUE) {
                rl_ok <- names(cuma[1:which(cuma > 99)[1]])
                # rl_ok<-as.character(seq(sort(as.numeric(rl_ok),decreasing =
                # FALSE)[1],sort(as.numeric(rl_ok),decreasing = TRUE)[1]))
                ps_comp <- ps_comp[names(ps_comp) %in% rl_ok]
            }
            
            signal_ps <- list()
            signal_ps_nt <- list()
            
            all_ps <- unlist(ps_comp)
            ps_comp[["all"]] <- all_ps
            mp <- mapToTranscripts(all_ps, transcripts = GTF_annotation$exons_txs[as.vector(seqnames(tile_cds))])
            mp$score <- all_ps$score[mp$xHits]
            # put seqlevels for compartment only to speed up stuff
            seqlevels(mp) <- as.vector(seqnames(tile_cds))
            seqlengths(mp) <- sum(width(GTF_annotation$exons_txs[as.vector(seqnames(tile_cds))]))
            covtx <- coverage(mp, weight = mp$score)
            ok_txs <- names(covtx[elementNROWS(runValue(covtx)) > 1])
            if (length(ok_txs) < 5) {
                ps_signals_tiles_all[[comp]] <- c()
                ps_signals_win_all[[comp]] <- c()
                next
            }
            # keep all txs if small n of txs
            
            if (length(ok_txs) >= 5 & length(covtx) < 250) {
                ok_txs <- names(covtx)
            }
            
            if (length(ok_txs) > 5000) {
                cnts_txss_agg <- sort(sum(covtx), decreasing = TRUE)
                ok_txs <- names(cnts_txss_agg)[1:5000]
            }
            if (fast_mode[bammo] == TRUE) {
                if (length(ok_txs) > 500) {
                  cnts_txss_agg <- sort(sum(covtx), decreasing = TRUE)
                  ok_txs <- names(cnts_txss_agg)[1:500]
                }
            }
            
            
            
            fivs <- GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords), ranges = IRanges(start = 1, 
                end = start(GTF_annotation$cds_txs_coords)), strand = "*")
            fivs <- fivs[as.character(seqnames(fivs)) %in% ok_txs]
            fivs_gen <- unlist(pmapFromTranscripts(fivs, transcripts = GTF_annotation$exons_txs[as.character(seqnames(fivs))], 
                ignore.strand = FALSE))
            fivs_gen <- fivs_gen[fivs_gen$hit]
            fivs_gen <- split(fivs_gen, names(fivs_gen))
            
            
            threes <- GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords), 
                ranges = IRanges(start = end(GTF_annotation$cds_txs_coords), end = GTF_annotation$cds_txs_coords$lentx), 
                strand = "*")
            threes <- threes[as.character(seqnames(threes)) %in% ok_txs]
            threes_gen <- unlist(pmapFromTranscripts(threes, transcripts = GTF_annotation$exons_txs[as.character(seqnames(threes))], 
                ignore.strand = FALSE))
            threes_gen <- threes_gen[threes_gen$hit]
            threes_gen <- split(threes_gen, names(threes_gen))
            
            cds_gen <- GTF_annotation$cds_txs[ok_txs]
            
            
            tile_cds <- tile_cds[as.character(seqnames(tile_cds)) %in% ok_txs]
            ok_txs <- unique(as.character(seqnames(tile_cds)))
            seqlevels(tile_cds) <- ok_txs
            list_covs <- list()
            list_covs[["all"]] <- covtx
            
            tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, end = start(tile_cds)))
            tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
            no5utr <- width(tile_5) < 51
            no3utr <- width(tile_3) < 51
            
            ex_annot <- GTF_annotation$exons_txs
            # if no utrs put bins to 1nt
            
            if (sum(no5utr) > 0) {
                tx_notok <- seqnames(tile_5)[no5utr]
                annot_notok <- ex_annot[tx_notok]
                annot_ok <- GRangesList(lapply(annot_notok, function(x) {
                  if (length(x) == 0) {
                    return(x)
                  }
                  x[1] <- resize(x[1], width = width(x[1]) + 51, fix = "end")
                  x
                }))
                ex_annot[names(annot_ok)] <- annot_ok
                
                seqlengths(tile_cds)[as.vector(tx_notok)] <- sum(width(annot_ok))
                tile_cds[no5utr] <- shift(tile_cds[no5utr], +51)
                tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, 
                  end = start(tile_cds)))
                tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                  end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                
                mp <- mapToTranscripts(all_ps, transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                seqlevels(mp) <- seqlevels(tile_cds)
                seqlengths(mp) <- seqlengths(tile_cds)
                mp$score <- all_ps$score[mp$xHits]
                
                covtx <- coverage(mp, weight = mp$score)
                list_covs[["all"]] <- covtx
            }
            
            if (sum(no3utr) > 0) {
                tx_notok <- seqnames(tile_3)[no3utr]
                annot_notok <- ex_annot[tx_notok]
                annot_ok <- GRangesList(lapply(annot_notok, function(x) {
                  if (length(x) == 0) {
                    return(x)
                  }
                  x[length(x)] <- resize(x[length(x)], width = width(x[length(x)]) + 
                    51, fix = "start")
                  x
                }))
                ex_annot[names(annot_ok)] <- annot_ok
                
                seqlengths(tile_cds)[as.vector(tx_notok)] <- sum(width(annot_ok))
                
                tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, 
                  end = start(tile_cds)))
                tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                  end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                
                mp <- mapToTranscripts(all_ps, transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                seqlevels(mp) <- seqlevels(tile_cds)
                seqlengths(mp) <- seqlengths(tile_cds)
                mp$score <- all_ps$score[mp$xHits]
                
                covtx <- coverage(mp, weight = mp$score)
                list_covs[["all"]] <- covtx
                
            }
            
            ps_tiles <- DataFrameList()
            ps_win <- DataFrameList()
            
            for (len in c("all", names(ps_comp))) {
                
                if ((sum(no3utr) + sum(no5utr)) == 0) {
                  mp <- mapToTranscripts(ps_comp[[len]], transcripts = fivs_gen)
                  mp$score <- ps_comp[[len]]$score[mp$xHits]
                  seqlevels(mp) <- names(fivs_gen)
                  seqlengths(mp) <- sum(width(fivs_gen))
                  cov_5 <- coverage(mp, weight = mp$score)
                  
                  
                  mp <- mapToTranscripts(ps_comp[[len]], transcripts = threes_gen)
                  mp$score <- ps_comp[[len]]$score[mp$xHits]
                  seqlevels(mp) <- names(threes_gen)
                  seqlengths(mp) <- sum(width(threes_gen))
                  cov_3 <- coverage(mp, weight = mp$score)
                  
                  mp <- mapToTranscripts(ps_comp[[len]], transcripts = cds_gen)
                  mp$score <- ps_comp[[len]]$score[mp$xHits]
                  seqlevels(mp) <- names(cds_gen)
                  seqlengths(mp) <- sum(width(cds_gen))
                  cov_cds <- coverage(mp, weight = mp$score)
                }
                
                if ((sum(no3utr) + sum(no5utr)) > 0) {
                  if (len != "all") {
                    mp <- mapToTranscripts(ps_comp[[len]], transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                    mp$score <- ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp) <- seqlevels(tile_cds)
                    seqlengths(mp) <- seqlengths(tile_cds)
                    covtx <- coverage(mp, weight = mp$score)
                    covtx <- covtx[ok_txs]
                  }
                  
                  cov_5 <- covtx[tile_5]
                  cov_3 <- covtx[tile_3]
                  cov_cds <- covtx[tile_cds]
                  
                  
                }
                
                ps_tiles_5 <- DataFrame(t(sapply(cov_5, function(x) {
                  clos <- 50 * (round(length(x)/50, digits = 0) + 1)
                  idx <- as.integer(seq(1, length(x), length.out = clos))
                  colMeans(matrix(x[idx], ncol = 50))
                })))
                ps_win_5 <- DataFrame(t(sapply(cov_5, function(x) {
                  as.vector(x)[c(1:25, (length(x) - 25):(length(x) - 1))]
                })))
                
                ps_tiles_3 <- DataFrame(t(sapply(cov_3, function(x) {
                  clos <- 50 * (round(length(x)/50, digits = 0) + 1)
                  idx <- as.integer(seq(1, length(x), length.out = clos))
                  colMeans(matrix(x[idx], ncol = 50))
                })))
                ps_win_3 <- DataFrame(t(sapply(cov_3, function(x) {
                  as.vector(x)[c(2:26, (length(x) - 24):length(x))]
                })))
                
                
                ps_tiles_cds <- DataFrame(t(sapply(cov_cds, function(x) {
                  clos <- 100 * (round(length(x)/100, digits = 0) + 1)
                  idx <- as.integer(seq(1, length(x), length.out = clos))
                  colMeans(matrix(x[idx], ncol = 100))
                })))
                ps_win_cds <- DataFrame(t(sapply(cov_cds, function(x) {
                  rnd <- as.integer(length(x)/2)%%3
                  mid <- (as.integer(length(x)/2) - rnd)
                  as.vector(x)[c(1:33, (mid - 17):(mid + 15), (length(x) - 32):length(x))]
                })))
                tls <- cbind(ps_tiles_5, ps_tiles_cds, ps_tiles_3)
                colnames(tls) <- c(paste("5_UTR", 1:length(ps_tiles_5[1, ]), sep = "_"), 
                  paste("CDS", 1:length(ps_tiles_cds[1, ]), sep = "_"), paste("3_UTR", 
                    1:length(ps_tiles_3[1, ]), sep = "_"))
                nts <- cbind(ps_win_5, ps_win_cds, ps_win_3)
                colnames(nts) <- c(paste("5_UTR", 1:length(ps_win_5[1, ]), sep = "_"), 
                  paste("CDS", 1:length(ps_win_cds[1, ]), sep = "_"), paste("3_UTR", 
                    1:length(ps_win_3[1, ]), sep = "_"))
                # DataFrameList?
                ps_tiles[[len]] <- tls
                ps_win[[len]] <- nts
                
            }
            
            
            ps_signals_tiles_all[[comp]] <- ps_tiles
            ps_signals_win_all[[comp]] <- ps_win
        }
        
        profiles_fivepr <- list(ps_signals_tiles_all, ps_signals_win_all)
        names(profiles_fivepr) <- c("five_prime_bins", "five_prime_subcodon")
        save(profiles_fivepr, file = paste(dira, "profiles_fivepr", sep = "/"))
        
        message(paste("Building aggregate 5' profiles --- Done!", date()))
        
        message(paste("Calculating P-sites offsets ...", date()))
        
        
        # calculate cutoffs for each read length
        
        list_cutoff <- lapply(profiles_fivepr[["five_prime_subcodon"]], function(x) {
            list_res <- list()
            for (i in names(x)) {
                if (i == "all") {
                  lnmax = 100
                } else {
                  lnmax = as.numeric(i)
                }
                list_res[[i]] <- calc_cutoffs_from_profiles(x[[i]], length_max = lnmax)
                
            }
            list_res
        })
        
        summary_res <- lapply(list_cutoff, function(x) t(sapply(x, function(y) {
            res <- c(y$final_cutoff, y$frames_res[2])
            names(res) <- c("cutoff", "frame_preference")
            res
        })))
        
        
        # choose readlengths
        
        res_rls <- choose_readlengths(summary_res, choice = readlength_choice_method[bammo], 
            nt_signals = profiles_fivepr[["five_prime_subcodon"]])
        
        message(paste("Calculating P-sites offsets --- Done!", date()))
        selection_cutoffs <- list(res_rls, summary_res, list_cutoff)
        names(selection_cutoffs) <- c("results_choice", "results_cutoffs", "analysis_frame_cutoff")
        save(selection_cutoffs, file = paste(dira, "selection_cutoffs", sep = "/"))
        
        
        
        # extract P-sites positions given selection above
        
        param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE), 
            what = c("mapq"), tag = "MD")
        seqllll <- seqlevels(unlist(read_stats$reads_pos1))
        seqleee <- seqlengths(unlist(read_stats$reads_pos1))
        
        # if we've provided our own choice of offsets, we should repeat it for each
        # compartment, should we need to
        if (!is.null(offsets_df)) {
            if (is.data.frame(offsets_df)) {
                res_rls$nucl$final_choice <- offsets_df %>% mutate(read_length = as.character(read_length)) %>% 
                  DataFrame
            } else {
                stopifnot(names(res_rls) %in% names(offsets_df))
                for (comp in names(res_rls)) res_rls[[comp]]$final_choice <- offsets_df[[comp]]
            }
        }
        
        rl_cutoffs_comp <- lapply(res_rls, function(x) {
            x$final_choice
        })
        
        # rescue all read lengths
        if (rescue_all_rls[bammo] == TRUE) {
            all_rlsss <- as.numeric(names(read_stats$reads_pos1))
            for (rilo in c("nucl", circs)) {
                ralo <- rl_cutoffs_comp[[rilo]]
                if (is.null(ralo)) {
                  ralo <- DataFrame()
                }
                rlno2 <- all_rlsss[!all_rlsss %in% ralo$read_length]
                ralo2 <- DataFrame(read_length = rlno2, cutoff = 12)
                rl_cutoffs_comp[[rilo]] <- rbind(ralo, ralo2)
            }
        }
        
        
        # what to do with chunks: x present chunk, y old chunks (cumulative)
        reduc <- function(x, y) {
            # the P-sites of the chunks are summed after the last chunk, all at once (sum_chunks()): the
            # merge of each chunk with the sums of the chunks before it took longer with each chunk
            chunks_ps <- x[["chunks_ps"]]
            if (is.null(chunks_ps)) {
                chunks_ps <- list(x[c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm")])
            }
            chunks_ps <- c(chunks_ps, list(y[c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm")]))
            
            covall_plus <- x$coverage_all_plus + y$coverage_all_plus
            covall_min <- x$coverage_all_minus + y$coverage_all_minus
            
            covuni_plus <- x$coverage_uniq_plus + y$coverage_uniq_plus
            covuni_min <- x$coverage_uniq_minus + y$coverage_uniq_minus
            
            
            rang_jun <- x$junctions
            rang_jun$reads <- rang_jun$reads + y$junctions$reads
            rang_jun$unique_reads <- rang_jun$unique_reads + y$junctions$unique_reads
            
            list_res <- list(x[["P_sites_all"]], x[["P_sites_uniq"]], x[["P_sites_uniq_mm"]], rang_jun, 
                covall_plus, covall_min, covuni_plus, covuni_min, chunks_ps)
            names(list_res) <- c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm", 
                "junctions", "coverage_all_plus", "coverage_all_minus", "coverage_uniq_plus", 
                "coverage_uniq_minus", "chunks_ps")
            
            
            return(list_res)
        }
        
        # what to do with each chunk (read as alignment file)
        
        yiel <- function(x) {
            readGAlignments(x, param = param)
        }
        
        # operations on the chunk (here count reads and whatnot)
        
        mapp <- function(x) {
            # the same strands as the offsets: mapp() of the counting pass inverts them too
            if (strandedness == "inverse") {
                x <- invertStrand(x)
            }
            psites_of_chunk(x, rl_cutoffs_comp, circs, seqllll, seqleee, GTF_annotation$junctions, 
                with_coverage = TRUE)
        }
        message(paste("Calculating P-sites positions and junctions ...", date()))
        
        P_sites_stats <- reduce_by_yield(X = opts, YIELD = yiel, MAP = mapp, REDUCE = reduc,
            n_cores = n_cores)
        # two or more chunks: sum the P-sites of all chunks
        if (!is.null(P_sites_stats[["chunks_ps"]])) {
            for (ps_set in c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm")) {
                chunks <- lapply(P_sites_stats$chunks_ps, function(x) x[[ps_set]])
                P_sites_stats[[ps_set]] <- sum_chunks(chunks, seqllll, seqleee)
            }
            P_sites_stats$chunks_ps <- NULL
        }
        
        save(P_sites_stats, file = paste(dira, "P_sites_stats", sep = "/"))
        
        message(paste("Calculating P-sites positions and junctions --- Done!", date()))
        
        ps <- P_sites_stats$P_sites_all
        
        message(paste("Building aggregate P-sites profiles ...", date()))
        
        # here do meta and other plots on P-sites positions ()
        ps_signals_win_all <- list()
        ps_signals_tiles_all <- list()
        codons_win_all <- list()
        ps_codons_ratio <- list()
        ps_codons_counts <- list()
        
        as_codons_ratio <- list()
        as_codons_counts <- list()
        
        es_codons_ratio <- list()
        es_codons_counts <- list()
        
        if (length(ps) == 0) {
            warning("Not enough signal or low frame preference, skipped P-sites & codon occupancy calculation. Set readlength_choice_method = 'all' to ignore this warning and get P-sites positions in a new run")
        }
        if (length(ps) > 0) {
            
            for (comp in c("nucl", circs)) {
                tile_cds <- list_txs_ok[[comp]]
                # put positive?
                strand(tile_cds) <- "+"
                tile_cds <- tile_cds[width(tile_cds) > 100]
                checkgen <- TRUE
                if (comp == "nucl") {
                  if ((sum((start(tile_cds) > 50))/length(start(tile_cds))) > 0.33) {
                    tile_cds <- tile_cds[start(tile_cds) > 50]
                    checkgen <- FALSE
                  }
                  if (sum(seqlengths(tile_cds)[as.vector(seqnames(tile_cds))] - end(tile_cds) > 
                    100)/length(tile_cds) > 0.33) {
                    tile_cds <- tile_cds[seqlengths(tile_cds)[as.vector(seqnames(tile_cds))] - 
                      end(tile_cds) > 100]
                    checkgen <- FALSE
                  }
                  ps_comp <- GRangesList(lapply(ps, function(x) {
                    sort(x[!seqnames(x) %in% circs])
                  }))
                } else {
                  ps_comp <- GRangesList(lapply(ps, function(x) {
                    sort(x[seqnames(x) == comp])
                  }))
                }
                
                
                signal_ps <- list()
                signal_ps_nt <- list()
                
                all_ps <- unlist(ps_comp)
                ps_comp[["all"]] <- all_ps
                mp <- mapToTranscripts(all_ps, transcripts = GTF_annotation$exons_txs[as.vector(seqnames(tile_cds))])
                mp$score <- all_ps$score[mp$xHits]
                # put seqlevels for compartment only to speed up stuff
                seqlevels(mp) <- names(GTF_annotation$exons_txs)
                seqlengths(mp) <- sum(width(GTF_annotation$exons_txs))
                covtx <- coverage(mp, weight = mp$score)
                ok_txs <- names(covtx[elementNROWS(runValue(covtx)) > 1])
                if (length(ok_txs) < 5) {
                  ps_signals_tiles_all[[comp]] <- c()
                  ps_signals_win_all[[comp]] <- c()
                  next
                }
                # keep all txs if small n of txs
                
                if (length(ok_txs) >= 5 & length(covtx) < 250) {
                  # covtx has all transcripts: keep the selected ones, as the 5' profiles do
                  ok_txs <- names(covtx)[names(covtx) %in% as.vector(seqnames(tile_cds))]
                }
                
                if (length(ok_txs) > 5000) {
                  cnts_txss_agg <- sort(sum(covtx), decreasing = TRUE)
                  ok_txs <- names(cnts_txss_agg)[seq_len(5000)]
                }
                
                if (fast_mode[bammo] == TRUE) {
                  if (length(ok_txs) > 500) {
                    cnts_txss_agg <- sort(sum(covtx), decreasing = TRUE)
                    ok_txs <- names(cnts_txss_agg)[seq_len(500)]
                  }
                }
                
                
                fivs <- GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords), 
                  ranges = IRanges(start = 1, end = start(GTF_annotation$cds_txs_coords)), 
                  strand = "*")
                fivs <- fivs[as.character(seqnames(fivs)) %in% ok_txs]
                fivs_gen <- unlist(pmapFromTranscripts(fivs, transcripts = GTF_annotation$exons_txs[as.character(seqnames(fivs))], 
                  ignore.strand = FALSE))
                fivs_gen <- fivs_gen[fivs_gen$hit]
                fivs_gen <- split(fivs_gen, names(fivs_gen))
                
                
                threes <- GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords), 
                  ranges = IRanges(start = end(GTF_annotation$cds_txs_coords), end = GTF_annotation$cds_txs_coords$lentx), 
                  strand = "*")
                threes <- threes[as.character(seqnames(threes)) %in% ok_txs]
                threes_gen <- unlist(pmapFromTranscripts(threes, transcripts = GTF_annotation$exons_txs[as.character(seqnames(threes))], 
                  ignore.strand = FALSE))
                threes_gen <- threes_gen[threes_gen$hit]
                threes_gen <- split(threes_gen, names(threes_gen))
                
                cds_gen <- GTF_annotation$cds_txs[ok_txs]
                
                
                tile_cds <- tile_cds[as.character(seqnames(tile_cds)) %in% ok_txs]
                ok_txs <- unique(as.character(seqnames(tile_cds)))
                seqlevels(tile_cds) <- ok_txs
                list_covs <- list()
                list_covs[["all"]] <- covtx
                
                tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, 
                  end = start(tile_cds)))
                tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                  end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                no5utr <- width(tile_5) < 51
                no3utr <- width(tile_3) < 51
                
                ex_annot <- GTF_annotation$exons_txs
                
                if (sum(no5utr) > 0) {
                  tx_notok <- seqnames(tile_5)[no5utr]
                  annot_notok <- ex_annot[tx_notok]
                  annot_ok <- GRangesList(lapply(annot_notok, function(x) {
                    if (length(x) == 0) {
                      return(x)
                    }
                    x[1] <- resize(x[1], width = width(x[1]) + 51, fix = "end")
                    x
                  }))
                  ex_annot[names(annot_ok)] <- annot_ok
                  
                  seqlengths(tile_cds)[as.vector(tx_notok)] <- sum(width(annot_ok))
                  tile_cds[no5utr] <- shift(tile_cds[no5utr], +51)
                  tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, 
                    end = start(tile_cds)))
                  tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                    end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                  
                  mp <- mapToTranscripts(all_ps, transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                  seqlevels(mp) <- seqlevels(tile_cds)
                  seqlengths(mp) <- seqlengths(tile_cds)
                  mp$score <- all_ps$score[mp$xHits]
                  
                  covtx <- coverage(mp, weight = mp$score)
                  list_covs[["all"]] <- covtx
                }
                
                if (sum(no3utr) > 0) {
                  tx_notok <- seqnames(tile_3)[no3utr]
                  annot_notok <- ex_annot[tx_notok]
                  annot_ok <- GRangesList(lapply(annot_notok, function(x) {
                    if (length(x) == 0) {
                      return(x)
                    }
                    x[length(x)] <- resize(x[length(x)], width = width(x[length(x)]) + 
                      51, fix = "start")
                    x
                  }))
                  ex_annot[names(annot_ok)] <- annot_ok
                  
                  seqlengths(tile_cds)[as.vector(tx_notok)] <- sum(width(annot_ok))
                  
                  tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, 
                    end = start(tile_cds)))
                  tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                    end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                  
                  mp <- mapToTranscripts(all_ps, transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                  seqlevels(mp) <- seqlevels(tile_cds)
                  seqlengths(mp) <- seqlengths(tile_cds)
                  mp$score <- all_ps$score[mp$xHits]
                  
                  covtx <- coverage(mp, weight = mp$score)
                  list_covs[["all"]] <- covtx
                  
                }
                
                
                ps_tiles <- DataFrameList()
                ps_win <- DataFrameList()
                cod_win <- DataFrameList()
                ps_cod <- DataFrameList()
                ps_cod_rat <- DataFrameList()
                
                as_cod <- DataFrameList()
                as_cod_rat <- DataFrameList()
                es_cod <- DataFrameList()
                es_cod_rat <- DataFrameList()
                txs_seqq <- NULL
                
                # ps_comp has "all" at the end: do it once, first
                for (len in c("all", setdiff(names(ps_comp), "all"))) {
                  
                  if ((sum(no3utr) + sum(no5utr)) == 0) {
                    mp <- mapToTranscripts(ps_comp[[len]], transcripts = fivs_gen)
                    mp$score <- ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp) <- names(fivs_gen)
                    seqlengths(mp) <- sum(width(fivs_gen))
                    cov_5 <- coverage(mp, weight = mp$score)
                    
                    
                    mp <- mapToTranscripts(ps_comp[[len]], transcripts = threes_gen)
                    mp$score <- ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp) <- names(threes_gen)
                    seqlengths(mp) <- sum(width(threes_gen))
                    cov_3 <- coverage(mp, weight = mp$score)
                    
                    mp <- mapToTranscripts(ps_comp[[len]], transcripts = cds_gen)
                    mp$score <- ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp) <- names(cds_gen)
                    seqlengths(mp) <- sum(width(cds_gen))
                    strand(mp) <- "+"
                    cov_cds <- coverage(mp, weight = mp$score)
                    cov_cds_a <- suppressWarnings(coverage(shift(mp, shift = 3), 
                      weight = mp$score))
                    cov_cds_e <- suppressWarnings(coverage(shift(mp, shift = -3), 
                      weight = mp$score))
                  }
                  
                  if ((sum(no3utr) + sum(no5utr)) > 0) {
                    
                    mp <- mapToTranscripts(ps_comp[[len]], transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                    mp$score <- ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp) <- seqlevels(tile_cds)
                    seqlengths(mp) <- seqlengths(tile_cds)
                    strand(mp) <- "+"
                    covtx <- coverage(mp, weight = mp$score)
                    covtx <- covtx[ok_txs]
                    cov_5 <- covtx[tile_5]
                    cov_3 <- covtx[tile_3]
                    cov_cds <- covtx[tile_cds]
                    cov_cds_a <- suppressWarnings(coverage(shift(mp, shift = 3), 
                      weight = mp$score))[tile_cds]
                    cov_cds_e <- suppressWarnings(coverage(shift(mp, shift = -3), 
                      weight = mp$score))[tile_cds]
                    
                  }
                  
                  ps_tiles_5 <- DataFrame(t(sapply(cov_5, function(x) {
                    clos <- 50 * (round(length(x)/50, digits = 0) + 1)
                    idx <- as.integer(seq(1, length(x), length.out = clos))
                    colMeans(matrix(x[idx], ncol = 50))
                  })))
                  ps_win_5 <- DataFrame(t(sapply(cov_5, function(x) {
                    as.vector(x)[c(seq_len(25), (length(x) - 25):(length(x) - 1))]
                  })))
                  
                  
                  ps_tiles_3 <- DataFrame(t(sapply(cov_3, function(x) {
                    clos <- 50 * (round(length(x)/50, digits = 0) + 1)
                    idx <- as.integer(seq(1, length(x), length.out = clos))
                    colMeans(matrix(x[idx], ncol = 50))
                  })))
                  ps_win_3 <- DataFrame(t(sapply(cov_3, function(x) {
                    as.vector(x)[c(2:26, (length(x) - 24):length(x))]
                  })))
                  
                  
                  ps_tiles_cds <- DataFrame(t(sapply(cov_cds, function(x) {
                    clos <- 100 * (round(length(x)/100, digits = 0) + 1)
                    idx <- as.integer(seq(1, length(x), length.out = clos))
                    colMeans(matrix(x[idx], ncol = 100))
                  })))
                  ps_win_cds <- DataFrame(t(sapply(cov_cds, function(x) {
                    rnd <- as.integer(length(x)/2)%%3
                    mid <- (as.integer(length(x)/2) - rnd)
                    as.vector(x)[c(seq_len(33), (mid - 17):(mid + 15), (length(x) - 
                      32):length(x))]
                  })))
                  
                  as_win_cds <- DataFrame(t(sapply(cov_cds_a, function(x) {
                    rnd <- as.integer(length(x)/2)%%3
                    mid <- (as.integer(length(x)/2) - rnd)
                    as.vector(x)[c(seq_len(33), (mid - 17):(mid + 15), (length(x) - 
                      32):length(x))]
                  })))
                  
                  es_win_cds <- DataFrame(t(sapply(cov_cds_e, function(x) {
                    rnd <- as.integer(length(x)/2)%%3
                    mid <- (as.integer(length(x)/2) - rnd)
                    as.vector(x)[c(seq_len(33), (mid - 17):(mid + 15), (length(x) - 
                      32):length(x))]
                  })))
                  
                  # select txs, output codon usage (each read length has the same transcripts: extract
                  # their sequences once)
                  if (is.null(txs_seqq) || !identical(names(txs_seqq), names(cov_cds))) {
                    txs_seqq <- extractTranscriptSeqs(x = genome_seq, transcripts = GTF_annotation$cds_txs[names(cov_cds)])
                  }
                  
                  gco <- as.character(names(getGeneticCode("1")))
                  names(gco) <- as.character(getGeneticCode("1"))
                  
                  # Get genetic code for compartment
                  
                  gen_cods <- GTF_annotation$genetic_codes
                  circ_cods <- which(comp == rownames(gen_cods))
                  if (length(circ_cods) > 0) {
                    
                    
                    gco <- as.character(names(getGeneticCode(gen_cods$genetic_code[circ_cods])))
                    names(gco) <- as.character(getGeneticCode(gen_cods$genetic_code[circ_cods]))
                    
                  }
                  
                  rnd <- as.integer(width(txs_seqq)/2)%%3
                  st_pos <- rep(1, length(txs_seqq))
                  mid_pos <- (as.integer(width(txs_seqq)/2) - rnd) - 17
                  end_pos <- width(txs_seqq) - 32
                  
                  iters <- seq(0, 10)
                  cod_counts <- c()
                  psit_counts <- c()
                  asit_counts <- c()
                  esit_counts <- c()
                  # Calculate codon occurrence and occupancy for different CDS sections OUTPUT AA?
                  for (i in iters) {
                    a <- narrow(txs_seqq, start = st_pos + 3 * i, end = st_pos + 
                      (2 + 3 * i))
                    coood <- as.character(a)
                    at <- table(a)
                    cod_cntt <- rep(0, length(gco))
                    names(cod_cntt) <- gco
                    cod_cntt[names(at)] <- as.numeric(at)
                    cod_counts <- cbind(cod_counts, cod_cntt)
                    pt <- rowSums(as.matrix(ps_win_cds[, (1 + 3 * i):(3 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    psit_counts <- cbind(psit_counts, ps_cntt)
                    
                    
                    pt <- rowSums(as.matrix(as_win_cds[, (1 + 3 * i):(3 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    asit_counts <- cbind(asit_counts, ps_cntt)
                    
                    
                    pt <- rowSums(as.matrix(es_win_cds[, (1 + 3 * i):(3 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    esit_counts <- cbind(esit_counts, ps_cntt)
                    
                  }
                  # mid
                  for (i in iters) {
                    a <- narrow(txs_seqq, start = mid_pos + 3 * i, end = mid_pos + 
                      (2 + 3 * i))
                    coood <- as.character(a)
                    at <- table(a)
                    cod_cntt <- rep(0, length(gco))
                    names(cod_cntt) <- gco
                    cod_cntt[names(at)] <- as.numeric(at)
                    cod_counts <- cbind(cod_counts, cod_cntt)
                    pt <- rowSums(as.matrix(ps_win_cds[, (34 + 3 * i):(36 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    psit_counts <- cbind(psit_counts, ps_cntt)
                    
                    pt <- rowSums(as.matrix(as_win_cds[, (34 + 3 * i):(36 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    asit_counts <- cbind(asit_counts, ps_cntt)
                    
                    pt <- rowSums(as.matrix(es_win_cds[, (34 + 3 * i):(36 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    esit_counts <- cbind(esit_counts, ps_cntt)
                  }
                  # end
                  for (i in iters) {
                    a <- narrow(txs_seqq, start = end_pos + 3 * i, end = end_pos + 
                      (2 + 3 * i))
                    coood <- as.character(a)
                    at <- table(a)
                    cod_cntt <- rep(0, length(gco))
                    names(cod_cntt) <- gco
                    cod_cntt[names(at)] <- as.numeric(at)
                    cod_counts <- cbind(cod_counts, cod_cntt)
                    pt <- rowSums(as.matrix(ps_win_cds[, (67 + 3 * i):(69 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    psit_counts <- cbind(psit_counts, ps_cntt)
                    
                    pt <- rowSums(as.matrix(as_win_cds[, (67 + 3 * i):(69 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    asit_counts <- cbind(asit_counts, ps_cntt)
                    
                    pt <- rowSums(as.matrix(es_win_cds[, (67 + 3 * i):(69 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    esit_counts <- cbind(esit_counts, ps_cntt)
                    
                  }
                  
                  colnames(psit_counts) <- c(paste("first", seq(1, 11), sep = "_"), 
                    paste("mid", seq(1, 11), sep = "_"), paste("last", seq(1, 11), 
                      sep = "_"))
                  colnames(esit_counts) <- c(paste("first", seq(1, 11), sep = "_"), 
                    paste("mid", seq(1, 11), sep = "_"), paste("last", seq(1, 11), 
                      sep = "_"))
                  colnames(asit_counts) <- c(paste("first", seq(1, 11), sep = "_"), 
                    paste("mid", seq(1, 11), sep = "_"), paste("last", seq(1, 11), 
                      sep = "_"))
                  colnames(cod_counts) <- c(paste("first", seq(1, 11), sep = "_"), 
                    paste("mid", seq(1, 11), sep = "_"), paste("last", seq(1, 11), 
                      sep = "_"))
                  ratio_psit_cod <- psit_counts/cod_counts
                  ratio_asit_cod <- asit_counts/cod_counts
                  ratio_esit_cod <- esit_counts/cod_counts
                  
                  psit_counts <- DataFrame(psit_counts)
                  asit_counts <- DataFrame(asit_counts)
                  esit_counts <- DataFrame(esit_counts)
                  cod_counts <- DataFrame(cod_counts)
                  ratio_psit_cod <- DataFrame(ratio_psit_cod)
                  ratio_asit_cod <- DataFrame(ratio_asit_cod)
                  ratio_esit_cod <- DataFrame(ratio_esit_cod)
                  
                  gcall <- paste(gco, names(gco), sep = ";")
                  
                  #remove 'Ns' our counts
                  non_n_rows <- grep(x=rownames(psit_counts),value=TRUE,pattern='N',invert=T)
                  psit_counts <- psit_counts[non_n_rows,]
                  asit_counts <- asit_counts[non_n_rows,]
                  esit_counts <- esit_counts[non_n_rows,]
                  ratio_psit_cod <- ratio_psit_cod[non_n_rows,]
                  ratio_asit_cod <- ratio_asit_cod[non_n_rows,]
                  ratio_esit_cod <- ratio_esit_cod[non_n_rows,]
                  cod_counts <- cod_counts[non_n_rows,]

                  rownames(psit_counts) <- gcall[match(rownames(psit_counts), gco)]
                  rownames(asit_counts) <- gcall[match(rownames(asit_counts), gco)]
                  rownames(esit_counts) <- gcall[match(rownames(esit_counts), gco)]
                  rownames(cod_counts) <- gcall[match(rownames(cod_counts), gco)]
                  
                  rownames(ratio_psit_cod) <- gcall[match(rownames(ratio_psit_cod), 
                    gco)]
                  rownames(ratio_asit_cod) <- gcall[match(rownames(ratio_asit_cod), 
                    gco)]
                  rownames(ratio_esit_cod) <- gcall[match(rownames(ratio_esit_cod), 
                    gco)]
                  
                  
                  cod_win[[len]] <- cod_counts
                  ps_cod[[len]] <- psit_counts
                  ps_cod_rat[[len]] <- ratio_psit_cod
                  
                  as_cod[[len]] <- asit_counts
                  as_cod_rat[[len]] <- ratio_asit_cod
                  
                  
                  es_cod[[len]] <- esit_counts
                  es_cod_rat[[len]] <- ratio_esit_cod
                  
                  
                  tls <- cbind(ps_tiles_5, ps_tiles_cds, ps_tiles_3)
                  colnames(tls) <- c(paste("5_UTR", seq_along(ps_tiles_5[1, ]), sep = "_"), 
                    paste("CDS", seq_along(ps_tiles_cds[1, ]), sep = "_"), paste("3_UTR", 
                      seq_along(ps_tiles_3[1, ]), sep = "_"))
                  nts <- cbind(ps_win_5, ps_win_cds, ps_win_3)
                  colnames(nts) <- c(paste("5_UTR", seq_along(ps_win_5[1, ]), sep = "_"), 
                    paste("CDS", seq_along(ps_win_cds[1, ]), sep = "_"), paste("3_UTR", 
                      seq_along(ps_win_3[1, ]), sep = "_"))
                  
                  
                  ps_tiles[[len]] <- tls
                  ps_win[[len]] <- nts
                  
                }
                
                codons_win_all[[comp]] <- cod_win
                ps_codons_ratio[[comp]] <- ps_cod_rat
                ps_codons_counts[[comp]] <- ps_cod
                
                as_codons_ratio[[comp]] <- as_cod_rat
                as_codons_counts[[comp]] <- as_cod
                
                
                es_codons_ratio[[comp]] <- es_cod_rat
                es_codons_counts[[comp]] <- es_cod
                
                ps_signals_tiles_all[[comp]] <- ps_tiles
                ps_signals_win_all[[comp]] <- ps_win
                
            }
        }
        
        profiles_P_sites <- list(ps_signals_tiles_all, ps_signals_win_all, codons_win_all, 
            ps_codons_counts, ps_codons_ratio, as_codons_counts, as_codons_ratio, 
            es_codons_counts, es_codons_ratio)
        
        names(profiles_P_sites) <- c("P_sites_bins", "P_sites_subcodon", "Codon_counts", 
            "P_sites_percodon", "P_sites_percodon_ratio", "A_sites_percodon", "A_sites_percodon_ratio", 
            "E_sites_percodon", "E_sites_percodon_ratio")
        save(profiles_P_sites, file = paste(dira, "profiles_P_sites", sep = "/"))
        # save(list =
        # c('ps_signals_tiles','ps_signals_tiles_all','ps_signals_win','ps_signals_win_all'),file
        # = 'ps_results_preprjan15')
        message(paste("Building aggregate P-sites profiles --- Done!", date()))
        
        message(paste("Exporting files ...", date()))
        
        merged_all_ps <- unlist(P_sites_stats$P_sites_all)
        if (length(merged_all_ps) > 0) {
            ps_strands <- sum_psites(merged_all_ps)
            
            merged_all_ps <- sort(c(ps_strands$plus, ps_strands$minus))
            
            export(ps_strands$plus, con = paste(dest_name, "_P_sites_plus.bedgraph", sep = ""))
            export(ps_strands$minus, con = paste(dest_name, "_P_sites_minus.bedgraph", sep = ""))
        }
        
        merged_uniq_ps <- unlist(P_sites_stats$P_sites_uniq)
        if (length(merged_uniq_ps) > 0) {
            
            ps_strands <- sum_psites(merged_uniq_ps)
            merged_uniq_ps <- sort(c(ps_strands$plus, ps_strands$minus))
            
            
            export(ps_strands$plus, con = paste(dest_name, "_P_sites_uniq_plus.bedgraph", 
                sep = ""))
            export(ps_strands$minus, con = paste(dest_name, "_P_sites_uniq_minus.bedgraph", 
                sep = ""))
        }
        
        
        merged_uniq_mm_ps <- unlist(P_sites_stats$P_sites_uniq_mm)
        if (length(merged_uniq_mm_ps) > 0) {
            
            ps_strands <- sum_psites(merged_uniq_mm_ps)
            merged_uniq_mm_ps <- sort(c(ps_strands$plus, ps_strands$minus))
        }
        
        if (!normalize_cov[bammo]) {
            export(P_sites_stats$coverage_all_plus, con = paste(dest_name, "_coverage_plus.bedgraph", 
                sep = ""))
            export(P_sites_stats$coverage_all_minus, con = paste(dest_name, "_coverage_minus.bedgraph", 
                sep = ""))
            export(P_sites_stats$coverage_uniq_plus, con = paste(dest_name, "_coverage_uniq_plus.bedgraph", 
                sep = ""))
            export(P_sites_stats$coverage_uniq_minus, con = paste(dest_name, "_coverage_uniq_minus.bedgraph", 
                sep = ""))
        }
        
        
        
        if (normalize_cov[bammo]) {
            bwpl <- P_sites_stats$coverage_all_plus
            bwmn <- P_sites_stats$coverage_all_minus
            
            correction <- (sum(as.numeric(unlist(runValue(bwpl)) * unlist(runLength(bwpl)))) + 
                sum(as.numeric(unlist(runValue(bwmn)) * unlist(runLength(bwmn)))))/1e+06
            okround <- sort(unique(unlist(runValue(bwpl))))[2]
            okround <- nchar(format(okround/correction, digits = 3, nsmall = 2)) - 
                2
            
            runValue(bwpl) <- round(as(runValue(bwpl), "NumericList")/correction, digits = okround)
            runValue(bwmn) <- round(as(runValue(bwmn), "NumericList")/correction, digits = okround)
            
            export(bwpl, con = paste(dest_name, "_coverage_plus.bedgraph", sep = ""))
            export(bwmn, con = paste(dest_name, "_coverage_minus.bedgraph", sep = ""))
            
            
            bwpl <- P_sites_stats$coverage_uniq_plus
            bwmn <- P_sites_stats$coverage_uniq_minus
            
            correction <- (sum(as.numeric(unlist(runValue(bwpl)) * unlist(runLength(bwpl)))) + 
                sum(as.numeric(unlist(runValue(bwmn)) * unlist(runLength(bwmn)))))/1e+06
            okround <- sort(unique(unlist(runValue(bwpl))))[2]
            okround <- nchar(format(okround/correction, digits = 3, nsmall = 2)) - 
                2
            
            runValue(bwpl) <- round(as(runValue(bwpl), "NumericList")/correction, digits = okround)
            runValue(bwmn) <- round(as(runValue(bwmn), "NumericList")/correction, digits = okround)
            
            export(bwpl, con = paste(dest_name, "_coverage_uniq_plus.bedgraph", sep = ""))
            export(bwmn, con = paste(dest_name, "_coverage_uniq_minus.bedgraph", 
                sep = ""))
            
        }
        
        
        juns <- P_sites_stats$junctions
        save(juns, file = paste(dest_name, "_junctions", sep = ""))
        
        
        for_ORFquant <- list(merged_all_ps, merged_uniq_ps, merged_uniq_mm_ps, juns)
        names(for_ORFquant) <- c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm", 
            "junctions")
        
        
        # Store top 50 occupied regions (possible contaminants)
        
        
        rds_st <- read_stats$reads_pos1
        rds_st <- sort(unlist(rds_st))
        
        regions <- list(reduce(unlist(GTF_annotation$cds_genes)), GTF_annotation$fiveutrs, 
            GTF_annotation$threeutrs, GTF_annotation$ncIsof, GTF_annotation$ncRNAs, 
            GTF_annotation$introns, GTF_annotation$intergenicRegions)
        names(regions) <- c("cds", "fiveutrs", "threeutrs", "ncIsof", "ncRNAs", "introns", 
            "intergenic")
        
        regions <- GRangesList(endoapply(regions, function(x) {
            mcols(x) <- NULL
            x
        }))
        
        rds_st <- rds_st[seqnames(rds_st) %in% unique(seqlevels(regions))]
        rds_st$len_adj <- as.numeric(names(rds_st))
        names(rds_st) <- NULL
        
        rds_st <- rds_st[order(rds_st$score, decreasing = TRUE)]
        rds_st$pct <- round(rds_st$score/(sum(rds_st$score)/100), digits = 4)
        top50 <- head(rds_st, 50)
        top50 <- resize(top50, width = top50$len_adj, fix = "start")
        top50$seq <- getSeq(genome_seq, top50)
        top50$region <- CharacterList("")
        
        ovg <- findOverlaps(top50, regions)
        ovg <- IntegerList(split(subjectHits(ovg), queryHits(ovg)))
        subsr <- as.numeric(names(ovg))
        a_reg <- CharacterList(lapply(ovg, FUN = function(x) {
            unique(names(regions)[x])
        }))
        top50$region[subsr] <- a_reg
        
        ovg <- findOverlaps(top50, GTF_annotation$genes)
        ovg <- IntegerList(split(subjectHits(ovg), queryHits(ovg)))
        subs <- as.numeric(names(ovg))
        a_nam <- CharacterList(lapply(ovg, FUN = function(x) {
            unique(names(GTF_annotation$genes)[x])
        }))
        a_biot <- CharacterList(lapply(a_nam, FUN = function(x) {
            GTF_annotation$trann$gene_biotype[match(x, GTF_annotation$trann$gene_id)]
        }))
        top50$gene <- CharacterList("")
        top50$gene_biotype <- CharacterList("")
        top50$gene[subs] <- a_nam
        top50$gene_biotype[subs] <- a_biot
        sequence_analysis <- top50
        
        # save as list_results
        
        save(sequence_analysis, file = paste(dira, "sequence_analysis", sep = "/"))
        
        
        res_all <- list(read_stats, profiles_fivepr, selection_cutoffs, P_sites_stats, 
            profiles_P_sites, sequence_analysis)
        
        names(res_all) <- c("read_stats", "profiles_fivepr", "selection_cutoffs", 
            "P_sites_stats", "profiles_P_sites", "sequence_analysis")
        if (length(merged_all_ps) > 0) {
            
            save(for_ORFquant, file = paste(dest_name, "_for_ORFquant", sep = ""))
        }
        finn <- lapply(res_all$selection_cutoffs$results_choice, function(x) {
            x$data
        })
        nope <- length(finn)
        
        namfinn <- rep(names(finn), elementNROWS(finn))
        finn <- do.call(rbind, args = finn)
        finn$comp <- namfinn
        
        pct_lib <- round(t(read_stats$rld/sum(read_stats$rld) * 100), digits = 4)
        rownames(pct_lib) <- gsub(rownames(pct_lib), pattern = "reads_", replacement = "")
        finn$pct_map <- 0
        for (i in unique(namfinn)) {
            finn$pct_map[finn$comp == i] <- pct_lib[match(finn$read_length[finn$comp == 
                i], rownames(pct_lib)), i]
        }
        if (nope == 0) {
            finn <- NULL
        }
        res_all$summary_P_sites <- finn
        if (nope == 0) {
            res_all$summary_P_sites <- c("")
        }
        if (write_tmp_files) {
            save(res_all, file = paste(dest_name, "_results_RiboseQC_all", sep = ""))
        }
        write.table(finn, file = paste(dest_name, "_P_sites_calcs", sep = ""), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
        # the offsets used for the P-sites, in the cutoff table format of ORFquant
        offs <- lapply(names(rl_cutoffs_comp), function(comp) {
            x <- rl_cutoffs_comp[[comp]]
            if (is.null(x) || nrow(x) == 0) {
                return(NULL)
            }
            data.frame(read_length = as.numeric(x$read_length), cutoff = as.numeric(x$cutoff),
                compartment = comp)
        })
        offs <- do.call(rbind, offs)
        if (is.null(offs)) {
            offs <- data.frame(read_length = numeric(), cutoff = numeric(), compartment = character())
        }
        write.table(offs, file = paste(dest_name, "_P_sites_offsets.tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        if (write_TSV_files) {
            write_tsv_files(res_all, finn, offs, dest_name, sample_names[bammo])
        }
        unlink(dira, recursive = TRUE)
        res_all$read_stats$reads_pos1 <- ""
        res_all$P_sites_stats <- ""
        
        nms <- names(res_all$profiles_fivepr$five_prime_bins)
        for (i in nms) {
            fvp <- res_all$profiles_fivepr$five_prime_bins[[i]]
            res_all$profiles_fivepr$five_prime_bins[[i]] <- lapply(fvp, function(x) {
                colSums(as.matrix(x))
            })
        }
        nms <- names(res_all$profiles_fivepr$five_prime_subcodon)
        for (i in nms) {
            fvp <- res_all$profiles_fivepr$five_prime_subcodon[[i]]
            res_all$profiles_fivepr$five_prime_subcodon[[i]] <- lapply(fvp, function(x) {
                colSums(as.matrix(x))
            })
        }
        
        nms <- names(res_all$profiles_P_sites$P_sites_bins)
        for (i in nms) {
            fvp <- res_all$profiles_P_sites$P_sites_bins[[i]]
            res_all$profiles_P_sites$P_sites_bins[[i]] <- lapply(fvp, function(x) {
                colSums(as.matrix(x))
            })
        }
        nms <- names(res_all$profiles_P_sites$P_sites_subcodon)
        for (i in nms) {
            fvp <- res_all$profiles_P_sites$P_sites_subcodon[[i]]
            res_all$profiles_P_sites$P_sites_subcodon[[i]] <- lapply(fvp, function(x) {
                colSums(as.matrix(x))
            })
        }
        resfile <- paste(dest_name, "_results_RiboseQC", sep = "")
        message(paste0('Writing results to ',resfile))
        save(res_all, file = resfile)
        rm(res_all, read_stats, profiles_fivepr, selection_cutoffs, P_sites_stats, 
            profiles_P_sites, sequence_analysis)
        gici <- gc()
        message(paste("Exporting files --- Done!", date(), "\n"))
        
        resfilelist <- c(resfilelist,resfile)
    }
    
      if (create_report) {
          message(paste("Creating html report in ", report_file, " ... ", date(), 
              sep = ""))
          
          create_html_report(input_files = resfilelist, 
              input_sample_names = sample_names, output_file = report_file, extended = extended_report)
          message(paste("Creating html report --- Done!", date(), "\n"))
          
          if (pdf_plots) {
              create_pdfs_from_rds_objects(output_rds_path = paste(report_file, "_plots/rds/", 
                  sep = ""))
          }
          
      }
    
    return(resfilelist)
    
}






# One range for each position of the ranges gr, with the score of its range. One chunk of a BAM file
# gives the 5' ends and the P-sites with one range for each position. The merge of two chunks adds
# their counts as coverage, so adjacent positions with the same count become one range, and the
# profiles lose its reads when it crosses the end of an exon or of the CDS (mapToTranscripts()
# gives no hit). With this, the merges give what one chunk gives.
single_positions <- function(gr) {
    wid <- width(gr)
    gr <- rep(gr, wid)
    ranges(gr) <- IRanges(start = start(gr) + sequence(wid) - 1L, width = 1L)
    gr
}

# The 5' ends or the P-sites of the chunks of a BAM file, summed: chunks is a list with a GRangesList of
# each chunk (a GRanges with the counts in "score" for each read length). The result has the read lengths
# in the order of the chunks, the chromosomes seql with the lengths seqle and one range for each position.
# Before, RiboseQC_analysis() merged each chunk with the sums of the chunks before it, and each merge
# took longer than the one before. One merge of all chunks gives the same sums.
sum_chunks <- function(chunks, seql, seqle) {
    all_ps <- GRangesList()
    for (rl in unique(unlist(lapply(chunks, names)))) {
        reads <- lapply(unname(chunks), function(x) {
            if (rl %in% names(x)) x[[rl]]
        })
        reads <- reads[lengths(reads) > 0]
        if (length(reads) > 0) {
            reads <- do.call(c, reads)
        } else {
            reads <- GRanges()
        }
        seqlevels(reads) <- seql
        seqlengths(reads) <- seqle
        pl <- reads[strand(reads) == "+"]
        mn <- reads[strand(reads) == "-"]
        if (length(pl) > 0) {
            covv_pl <- coverage(pl, weight = pl$score)
        } else {
            covv_pl <- coverage(pl)
        }
        covv_pl <- GRanges(covv_pl)
        covv_pl <- covv_pl[covv_pl$score > 0]
        if (length(mn) > 0) {
            covv_min <- coverage(mn, weight = mn$score)
        } else {
            covv_min <- coverage(mn)
        }
        covv_min <- GRanges(covv_min)
        covv_min <- covv_min[covv_min$score > 0]
        strand(covv_pl) <- "+"
        strand(covv_min) <- "-"
        all_ps[[rl]] <- single_positions(sort(c(covv_pl, covv_min)))
    }
    all_ps
}

# The tab-separated files of RiboseQC_analysis(write_TSV_files = TRUE) for one BAM file, with the prefix
# dest_name: tables of res_all, and a summary with one row for each compartment. finn is the table of
# P_sites_calcs (NULL when no read length has an offset), offs the table of P_sites_offsets.tsv
write_tsv_files <- function(res_all, finn, offs, dest_name, sample_name) {
    write_tsv <- function(x, suffix) {
        write.table(x, file = paste(dest_name, suffix, sep = ""), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
    }
    rs <- res_all$read_stats
    write_tsv(data.frame(compartment = rownames(rs$rld), rs$rld, check.names = FALSE), "_read_lengths.tsv")
    write_tsv(data.frame(region = rownames(rs$positions), rs$positions, check.names = FALSE), "_regions.tsv")
    write_tsv(data.frame(gene_id = rownames(rs$counts_all_genes), as.data.frame(rs$counts_all_genes),
        check.names = FALSE), "_genes.tsv")
    write_tsv(finn, "_P_sites_calcs.tsv")

    # the ranges, then their columns: sequences as text, lists with their values separated by commas
    top50 <- res_all$sequence_analysis
    top <- data.frame(seqnames = as.character(seqnames(top50)), start = start(top50), end = end(top50),
        width = width(top50), strand = as.character(strand(top50)))
    for (n in names(mcols(top50))) {
        v <- mcols(top50)[[n]]
        if (is(v, "XStringSet")) {
            v <- as.character(v)
        } else if (is(v, "List")) {
            v <- vapply(as.list(v), paste, "", collapse = ",")
        }
        top[[n]] <- v
    }
    write_tsv(top, "_top50_positions.tsv")

    # the reads of the read lengths with an offset give the P-sites; the frame preference is the mean of
    # these read lengths, weighted by their reads, without those that have none in P_sites_calcs
    fp_all <- numeric()
    if (!is.null(finn)) {
        fp_all <- finn$frame_preference
        names(fp_all) <- paste(finn$comp, finn$read_length)
    }
    summ <- lapply(rownames(rs$rld), function(comp) {
        rl_reads <- unlist(rs$rld[comp, ])
        names(rl_reads) <- sub("^reads_", "", names(rl_reads))
        reads <- sum(rl_reads)
        rls <- as.character(offs$read_length[offs$compartment == comp])
        ps_reads <- rl_reads[rls]
        ps_reads[is.na(ps_reads)] <- 0
        fp <- fp_all[paste(comp, rls)]
        ok <- !is.na(fp) & ps_reads > 0
        pct <- function(x) if (reads > 0) round(x/reads * 100, digits = 4) else NA
        data.frame(sample = sample_name, compartment = comp, reads = reads,
            pct_cds = pct(rs$positions["cds", paste("reads", comp, sep = "_")]), pct_P_sites = pct(sum(ps_reads)),
            frame_preference = if (any(ok)) round(sum(fp[ok] * ps_reads[ok])/sum(ps_reads[ok]), digits = 4) else NA)
    })
    write_tsv(do.call(rbind, summ), "_summary.tsv")
}

# Read the BamFile X in chunks with YIELD, run MAP on each chunk and combine the results with
# REDUCE, in the order of the chunks. This is the loop of GenomicFiles::reduceByYield(), which
# RiboseQC used before. With n_cores > 1, forked processes run MAP on the chunks, at most n_cores
# at a time, while the next ones are read (the loop of prepare_for_ORFquant(), from ORFquant).
reduce_by_yield <- function(X, YIELD, MAP, REDUCE, n_cores = 1) {
    if (!isOpen(X)) {
        open(X)
        on.exit(close(X))
    }
    result <- list()
    n_chunks <- 0
    jobs <- list()
    eof <- FALSE
    while (!eof || length(jobs) > 0) {
        if (!eof && length(jobs) < n_cores) {
            data <- YIELD(X)
            eof <- length(data) == 0L
            if (!eof)
                jobs[[length(jobs) + 1]] <- mc_job(MAP(data), n_cores)
            rm(data)
        }
        # the results of the oldest jobs that are done, so that their processes end; wait for the
        # oldest when n_cores are running or all chunks are read
        while (length(jobs) > 0) {
            res <- mc_value(jobs[[1]], wait = eof || length(jobs) >= n_cores)
            if (is.null(res))
                break
            jobs <- jobs[-1]
            n_chunks <- n_chunks + 1
            if (n_chunks == 1) {
                result <- res
            } else {
                result <- REDUCE(result, res)
            }
        }
    }
    result
}

# RiboseQC_analysis() does not use its argument genome_seq (a path or an FaFile): the genome sequence is
# annotation_genome, the FASTA file of the annotation. This warns when genome_seq is another file
check_genome_seq <- function(genome_seq, annotation_genome, annotation_file) {
    if (is.null(genome_seq))
        return(invisible(NULL))
    given <- if (is(genome_seq, "FaFile")) path(genome_seq) else genome_seq
    annot <- if (is(annotation_genome, "FaFile")) path(annotation_genome) else NA_character_
    if (!is.character(given) || length(given) != 1 || is.na(annot) ||
        normalizePath(given, mustWork = FALSE) != normalizePath(annot, mustWork = FALSE)) {
        warning("genome_seq (", if (is.character(given)) paste(given, collapse = ", ") else class(genome_seq)[1],
            ") is not used: RiboseQC_analysis() uses the genome sequence of the annotation ", annotation_file,
            if (!is.na(annot)) paste0(", the FASTA file ", annot), ". To use another FASTA file, make the annotation with it",
            call. = FALSE)
    }
}

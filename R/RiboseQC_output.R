#' Output files of RiboseQC_analysis() and the list res_all
#'
#' \code{\link{RiboseQC_analysis}} writes files for each BAM file, with its entry of \code{dest_names} as the prefix,
#' and the html report of all BAM files. The file \code{<dest_name>_results_RiboseQC_all} has all results, in the list
#' \code{res_all}. This page describes each file, each element of \code{res_all}, the files of the report and the file
#' of \code{\link{prepare_for_ORFquant}}. \code{\link{prepare_annotation_files}} describes the files of the annotation.
#'
#' @details
#' The words of this page:
#' \describe{
#' \item{compartment}{\code{"nucl"}, all chromosomes that are not circular (also those that are not in the
#'   annotation), and each circular chromosome of the annotation (\code{circ_chroms} of
#'   \code{\link{prepare_annotation_files}}) that is in the header of the BAM file, with its name (for example
#'   \code{"chrM"}). RiboseQC calculates the offsets and the profiles for each compartment.}
#' \item{reads}{The alignments of the BAM file, without the alignments flagged as duplicates, the secondary
#'   alignments and the alignments with an insertion or a deletion. Each alignment counts as one read.}
#' \item{read length}{The length of the aligned part of the read: soft-clipped bases are not counted.}
#' \item{uniquely mapping reads}{Reads with a mapping quality of more than 50. Their elements and files have
#'   \code{uniq} or \code{_unq} in the name.}
#' \item{5' end}{The first aligned nt of the read, on its strand. With \code{stranded = "inverse"}, RiboseQC first
#'   inverts the strand of each read.}
#' \item{offset (cutoff)}{The distance of the P-site from the 5' end, in nt of the aligned blocks of the read (so a
#'   spliced read gets its P-site in the next exon when the offset goes past the end of its first block). A read of
#'   a read length that has an offset in its compartment gives one P-site.}
#' \item{regions}{The 7 regions of the annotation (see \code{\link{prepare_annotation_files}}): \code{cds} (the
#'   CDS of \code{cds_genes}), \code{fiveutrs}, \code{threeutrs}, \code{ncIsof}, \code{ncRNAs}, \code{introns} and
#'   \code{intergenic} (\code{intergenicRegions}, without a strand).}
#' }
#'
#' @section Files of each BAM file:
#' The names are the entry of \code{dest_names} followed by the suffix, with no separator: \code{dest_names =
#' "results/root"} gives \code{results/root_results_RiboseQC}. The R files hold one object each, to read with
#' \code{load()}.
#' \describe{
#' \item{\code{_results_RiboseQC_all}}{The list \code{res_all}, described below. Only with \code{write_tmp_files =
#'   TRUE} (the default).}
#' \item{\code{_results_RiboseQC}}{A smaller \code{res_all}, for \code{\link{create_html_report}} (see "The smaller
#'   results file" below). Always written.}
#' \item{\code{_P_sites_calcs}}{The table \code{summary_P_sites}, tab-separated, with a header line. An empty line when
#'   no read length has an offset.}
#' \item{\code{_P_sites_offsets.tsv}}{The offsets that give the P-sites, tab-separated, with the columns
#'   \code{read_length}, \code{cutoff} and \code{compartment}. For each compartment, the read lengths that
#'   \code{readlength_choice_method} selects (or those of \code{offsets_df}), then the read lengths of
#'   \code{rescue_all_rls}. Only the header line when there is no offset. \code{\link{prepare_for_ORFquant}} reads
#'   this file.}
#' \item{\code{_for_ORFquant}}{The list \code{for_ORFquant}, for ORFquant: \code{P_sites_all}, \code{P_sites_uniq} and
#'   \code{P_sites_uniq_mm} (the P-sites of \code{P_sites_stats}, summed over the read lengths and the compartments)
#'   and \code{junctions} (the same as in \code{P_sites_stats}). Each P-site GRanges has the number of P-sites at each
#'   nt in \code{score}. Neighbouring positions with the same number are one range, so a range can be longer than
#'   1 nt, and the sum of \code{score} is less than the number of P-sites when a range is longer
#'   (\code{sum(score * width)} is the number). Only written when there are P-sites.}
#' \item{\code{_junctions}}{The GRanges \code{juns}: \code{junctions} of \code{P_sites_stats}.}
#' \item{\code{_P_sites_plus.bedgraph}, \code{_P_sites_minus.bedgraph}, \code{_P_sites_uniq_plus.bedgraph},
#'   \code{_P_sites_uniq_minus.bedgraph}}{The ranges of \code{P_sites_all} and \code{P_sites_uniq} of
#'   \code{for_ORFquant} on each strand, in the bedGraph format (start 0-based, without a track line). The two
#'   files of all reads are written when there are P-sites, and the two files of uniquely mapping reads when these
#'   have P-sites. A strand without P-sites gives an empty file. There is no file of \code{P_sites_uniq_mm}.}
#' \item{\code{_coverage_plus.bedgraph}, \code{_coverage_minus.bedgraph}, \code{_coverage_uniq_plus.bedgraph},
#'   \code{_coverage_uniq_minus.bedgraph}}{The coverage of the reads on each strand (\code{coverage_*} of
#'   \code{P_sites_stats}: the aligned blocks, without the introns and the soft-clipped bases), for all chromosomes
#'   of the BAM header, also the ranges with no coverage. With \code{normalize_cov = TRUE} (the default), the
#'   coverage is divided by the coverage of the two strands together (the number of aligned nt) and multiplied by
#'   1 million, and rounded. The files of the uniquely mapping reads use their own total. With \code{normalize_cov =
#'   FALSE}, the coverage is the number of reads.}
#' \item{\code{_summary.tsv}, \code{_read_lengths.tsv}, \code{_regions.tsv}, \code{_genes.tsv},
#'   \code{_P_sites_calcs.tsv}, \code{_top50_positions.tsv}}{The tab-separated files of \code{write_TSV_files = TRUE}
#'   (the default), see "Tab-separated files" below.}
#' }
#' During the analysis of a BAM file, RiboseQC saves parts of the results in the folder
#' \code{tmp_RiboseQC_<name>} (\code{<name>} is the last part of \code{dest_name}), in the folder of
#' \code{dest_name}. It deletes the folder before it saves \code{_results_RiboseQC}. The folder stays if the function
#' stops with an error.
#'
#' @section Tab-separated files:
#' With \code{write_TSV_files = TRUE}, the tables of \code{res_all} are also written as tab-separated files with a
#' header line, with the same names as in \code{res_all}:
#' \describe{
#' \item{\code{_summary.tsv}}{One row for each compartment, with the columns \code{sample} (of \code{sample_names}),
#'   \code{compartment}, \code{reads}, \code{pct_cds} (the percentage of the reads of the compartment on its
#'   \code{cds} region), \code{pct_P_sites} (the percentage of the reads of the compartment with a read length that
#'   has an offset in \code{_P_sites_offsets.tsv}: these reads give the P-sites) and \code{frame_preference} (the mean
#'   \code{frame_preference} of these read lengths, weighted by their reads). Read lengths without a
#'   \code{frame_preference} in \code{_P_sites_calcs}, from \code{offsets_df} or \code{rescue_all_rls}, are not in
#'   the mean; it is \code{NA} if no read length has one. The percentages are \code{NA} for a compartment without
#'   reads.}
#' \item{\code{_read_lengths.tsv}}{\code{read_stats$rld}: one row for each compartment, with the column
#'   \code{compartment} and the number of reads of each read length (\code{reads_<read length>}).}
#' \item{\code{_regions.tsv}}{\code{read_stats$positions}: one row for each region, with the column \code{region} and
#'   the number of reads on it in each compartment (\code{reads_<compartment>}). A read on two regions counts for
#'   both.}
#' \item{\code{_genes.tsv}}{\code{read_stats$counts_all_genes}: one row for each gene, with the columns
#'   \code{gene_id}, \code{reads}, \code{chromosome}, \code{gene_name}, \code{gene_biotype}, \code{RPKM} and
#'   \code{TPM}.}
#' \item{\code{_P_sites_calcs.tsv}}{The same table as \code{_P_sites_calcs} (\code{summary_P_sites}).}
#' \item{\code{_top50_positions.tsv}}{\code{sequence_analysis}: one row for each of its 50 ranges, with the columns
#'   \code{seqnames}, \code{start}, \code{end}, \code{width} and \code{strand} of the range, then \code{score},
#'   \code{len_adj}, \code{pct}, \code{seq}, \code{region}, \code{gene} and \code{gene_biotype}. More values are
#'   separated by commas.}
#' }
#'
#' @section res_all$read_stats:
#' The reads, by compartment, read length, region and gene.
#' \describe{
#' \item{\code{rld}}{A data.frame with one row for each compartment (the row names) and one column
#'   \code{reads_<read length>} for each read length from the shortest to the longest read of the BAM file (also the
#'   read lengths without reads): the number of reads.}
#' \item{\code{rld_unq}}{The same for the uniquely mapping reads.}
#' \item{\code{positions}}{A data.frame with one row for each region (the row names) and one column
#'   \code{reads_<compartment>} for each compartment: the number of reads that overlap the region by at least 1 nt
#'   (the aligned blocks of the read, not only its 5' end). A read on two regions counts for both, and a read on no
#'   region is not counted. With \code{stranded = FALSE}, the strand of the regions is ignored.}
#' \item{\code{positions_unq}}{The same for the uniquely mapping reads.}
#' \item{\code{reads_pos1}}{A GRangesList with one GRanges for each read length with reads (the names): the 5' ends of
#'   the reads of all compartments, 1 nt wide, with the number of reads in \code{score}. The sum of \code{score} of a
#'   read length is the sum of its column of \code{rld}.}
#' \item{\code{counts_cds_genes}}{A DataFrame with one row for each gene with a CDS (the gene ids are the row names),
#'   and the columns \code{reads} (the reads that overlap the CDS of the gene by at least 1 nt; a read counts for
#'   each gene that it overlaps; strand as in \code{positions}), \code{chromosome}, \code{gene_name} and
#'   \code{gene_biotype} (from the GTF file, \code{NA} when it has none), \code{RPKM} and \code{TPM}. The length of a
#'   gene is the width of its CDS (all its CDS together). \code{RPKM} is \code{reads} per 1,000 nt of length and per
#'   million of the sum of \code{reads} of the table; \code{TPM} is \code{reads} per 1,000 nt of length, scaled to
#'   a sum of 1 million. Both are rounded to 4 digits.}
#' \item{\code{counts_all_genes}}{The same for all genes, with their exons instead of their CDS.}
#' \item{\code{counts_cds_genes_unq}, \code{counts_all_genes_unq}}{The same for the uniquely mapping reads.}
#' \item{\code{reads_summary}}{A list with one DataFrame for each compartment, with one row for each region and the
#'   columns of \code{rld}: the reads of \code{positions}, by read length (strand as in \code{positions}).
#'   \code{read_subset} selects the read lengths with it.}
#' \item{\code{reads_summary_unq}}{The same for the uniquely mapping reads.}
#' }
#'
#' @section res_all$profiles_fivepr:
#' The 5' ends of the reads on the transcripts, for each compartment and read length.
#' \describe{
#' \item{\code{five_prime_bins}}{A list with one element for each compartment with profiles (the names). Each element
#'   is a DataFrameList with \code{"all"} and one DataFrame for each read length: the read lengths that
#'   \code{read_subset} selects (with \code{read_subset = FALSE}, all read lengths), from the shortest to the longest. \code{"all"} is
#'   the sum of these read lengths. A DataFrame has one row for each transcript (the transcript ids are the row
#'   names) and 200 columns, \code{5_UTR_1} to \code{5_UTR_50}, \code{CDS_1} to \code{CDS_100} and \code{3_UTR_1} to
#'   \code{3_UTR_50}. The 5' part (from the first nt of the transcript to the first nt of the start codon), the CDS
#'   and the 3' part (from the last nt of the CDS to the last nt of the transcript) are each cut into 50, 100 and 50
#'   bins of about the same length. Each value is the mean number of 5' ends for each nt of the bin.}
#' \item{\code{five_prime_subcodon}}{The same, with 199 columns: the number of 5' ends at each nt of 7 windows.
#'   \code{5_UTR_1} to \code{5_UTR_25}: the first 25 nt of the transcript. \code{5_UTR_26} to \code{5_UTR_50}: the
#'   25 nt before the start codon. \code{CDS_1} to \code{CDS_33}: the first 33 nt of the CDS (\code{CDS_1} is the
#'   first nt of the start codon). \code{CDS_34} to \code{CDS_66}: 33 nt in the middle of the CDS, from the first nt
#'   of a codon. \code{CDS_67} to \code{CDS_99}: the last 33 nt of the CDS (the last 3 are the stop codon when the
#'   CDS of the annotation includes it). \code{3_UTR_1} to \code{3_UTR_25}: the 25 nt after the CDS.
#'   \code{3_UTR_26} to \code{3_UTR_50}: the last 25 nt of the transcript.}
#' }
#' The transcripts: in \code{nucl}, the representative transcript of each gene (\code{reprentative_mostcommon} of
#' \code{cds_txs_coords}); in a circular compartment, all transcripts with a CDS. Only transcripts with a CDS longer
#' than 100 nt. In \code{nucl}, if more than a third of them have a 5'UTR of at least 50 nt, only these; and if more
#' than a third of them have a 3'UTR longer than 100 nt, only these. Of these transcripts, the rows are those with a
#' 5' end of the read lengths of the profiles. If fewer than 250 transcripts are left after the filters, the rows are
#' all of them, also those without 5' ends; if more than 5000 have 5' ends, the 5000 with the most (500 with
#' \code{fast_mode = TRUE}). A compartment
#' with fewer than 5 transcripts with 5' ends, or without reads on its CDS and UTRs, has no profiles. When a
#' transcript has a 5'UTR shorter than 50 nt, RiboseQC adds 51 nt of the genome before its first exon, and when it has
#' a 3'UTR shorter than 50 nt, 51 nt after its last exon, so its windows can include nt that are not in the
#' transcript.
#'
#' @section res_all$selection_cutoffs:
#' The offsets that RiboseQC calculates from \code{five_prime_subcodon} (see \code{\link{calc_cutoffs_from_profiles}}
#' and \code{\link{choose_readlengths}}). \code{offsets_df} and \code{rescue_all_rls} do not change them.
#' \describe{
#' \item{\code{results_choice}}{A list with one element for each compartment with at least one read length with an
#'   offset: \code{data}, a DataFrame with one row for each read length with an offset (the row names), by decreasing
#'   \code{frame_preference}, with the columns of \code{summary_P_sites} from \code{cutoff} to \code{max_coverage};
#'   and \code{final_choice}, the columns \code{read_length} and \code{cutoff} of the rows of \code{data} that
#'   \code{readlength_choice_method} selects.}
#' \item{\code{results_cutoffs}}{A list with one numeric matrix for each compartment with profiles: one row for
#'   \code{"all"} and for each read length with an offset (the row names), and the columns \code{cutoff} and
#'   \code{frame_preference}. A matrix with 1 row and no column when no read length has an offset.}
#' \item{\code{analysis_frame_cutoff}}{A list with one element for each compartment with profiles, with one element
#'   for \code{"all"} and for each read length with an offset: the value of \code{\link{calc_cutoffs_from_profiles}},
#'   a list of \code{final_cutoff} (the offset), \code{frames_res} (\code{max_frame_in_pctORF}, the percentage of
#'   transcripts with the most 5' ends in the frame that most transcripts have, and \code{mean_pct_max_frame}, the
#'   \code{frame_preference}), \code{cutoffs} (the offsets of the 8 methods, \code{NA} when a method cannot
#'   calculate one), \code{frames} (a matrix with one row for each transcript with 5' ends in the 99 nt of the CDS
#'   windows, and the fractions of these 5' ends in the first, the second and the third nt of the codons) and
#'   \code{profiles_start} (the 5' ends of each transcript from 24 nt before the start codon to its first nt,
#'   columns \code{"-24"} to \code{"0"}; for a read length shorter than 24, from the read length). An empty list
#'   when no read length has an offset.}
#' }
#'
#' @section res_all$P_sites_stats:
#' \describe{
#' \item{\code{P_sites_all}}{A GRangesList with one GRanges for each read length with an offset (the names, in the
#'   order of \code{_P_sites_offsets.tsv}): the P-sites of the reads of this read length in all compartments, each
#'   with the offset of its compartment, 1 nt wide, on the strand of the read, with the number of P-sites in
#'   \code{score}. The offsets are those of \code{_P_sites_offsets.tsv}, so \code{offsets_df} and
#'   \code{rescue_all_rls} change them.}
#' \item{\code{P_sites_uniq}}{The same for the uniquely mapping reads.}
#' \item{\code{P_sites_uniq_mm}}{The same for the uniquely mapping reads with a mismatch in the MD tag of the BAM file.
#'   Empty when the BAM file has no MD tags.}
#' \item{\code{junctions}}{A GRanges with one range for each intron of the annotation (\code{junctions} of the
#'   annotation), with the columns \code{type} (\code{"J"}), \code{tx_name} and \code{gene_id} (the transcripts and
#'   the genes of the intron) and \code{reads} and \code{unique_reads}: the number of reads (of all read lengths and
#'   compartments), and of uniquely mapping reads, with a gap (\code{N} in the CIGAR string) with the start, the end
#'   and the strand of the intron.}
#' \item{\code{coverage_all_plus}, \code{coverage_all_minus}, \code{coverage_uniq_plus},
#'   \code{coverage_uniq_minus}}{RleLists with the number of reads at each nt of each chromosome of the BAM header,
#'   on the plus and the minus strand, for all reads and for the uniquely mapping reads (the aligned blocks of the
#'   reads, without the introns and the soft-clipped bases).}
#' }
#'
#' @section res_all$profiles_P_sites:
#' The P-sites on the transcripts and on the codons, for each compartment with P-sites on at least 5 transcripts.
#' Each element is a list with one element for each such compartment, with \code{"all"} and one element for each read
#' length with an offset in any compartment, in the order of \code{_P_sites_offsets.tsv} (so a read length can have
#' only zeros in a compartment). \code{"all"} is the sum of the read lengths. The rows are the transcripts that
#' \code{profiles_fivepr} can use (see above) with at least one P-site (at most 5000, or 500 with
#' \code{fast_mode = TRUE}). With
#' \code{offsets_df} or \code{rescue_all_rls}, other P-sites can give other transcripts. Without P-sites, each element
#' is an empty list.
#' \describe{
#' \item{\code{P_sites_bins}, \code{P_sites_subcodon}}{The P-sites, with the rows and the columns of
#'   \code{five_prime_bins} and \code{five_prime_subcodon}.}
#' \item{\code{Codon_counts}}{A DataFrame with one row for each codon, named \code{<codon>;<amino acid>} (for example
#'   \code{"TTT;F"}, and \code{"TAA;*"} for a stop codon), and 33 columns for the 33 codons of the CDS windows of
#'   \code{P_sites_subcodon}: \code{first_1} to \code{first_11} (\code{first_1} is the start codon), \code{mid_1}
#'   to \code{mid_11} and \code{last_1} to \code{last_11} (\code{last_11} is the stop codon when the CDS of the
#'   annotation includes it). The value is the number of transcripts with this codon at this position, the same for
#'   all read lengths. The amino acids are those of the genetic code of the compartment: the standard code for
#'   \code{nucl}, and the code of the annotation (\code{genetic_codes}) for a circular chromosome.}
#' \item{\code{P_sites_percodon}}{The same layout: the number of P-sites on the 3 nt of the codon, at this position, in
#'   the transcripts that have this codon there.}
#' \item{\code{P_sites_percodon_ratio}}{\code{P_sites_percodon} divided by \code{Codon_counts}: the mean number of
#'   P-sites for each occurrence of the codon. \code{NaN} when the codon does not occur at this position.}
#' \item{\code{A_sites_percodon}, \code{A_sites_percodon_ratio}}{The same for the reads with the A site on the codon:
#'   the P-sites on the codon before it.}
#' \item{\code{E_sites_percodon}, \code{E_sites_percodon_ratio}}{The same for the reads with the E site on the codon:
#'   the P-sites on the codon after it.}
#' }
#'
#' @section res_all$sequence_analysis and res_all$summary_P_sites:
#' \describe{
#' \item{\code{sequence_analysis}}{A GRanges with the 50 combinations of a 5' end and a read length with the most
#'   reads, on the chromosomes of the annotation. Each range starts at the 5' end and has the width of the read
#'   length. The columns: \code{score} (the number of reads), \code{len_adj} (the read length), \code{pct} (the
#'   percentage of all reads on the chromosomes of the annotation, rounded to 4 digits), \code{seq} (the sequence of
#'   the genome in the range, on its strand), \code{region} (the regions that the range overlaps, \code{""} if none),
#'   \code{gene} and \code{gene_biotype} (the ids and the biotypes of the genes that the range overlaps on its strand,
#'   \code{""} if none).}
#' \item{\code{summary_P_sites}}{A DataFrame with one row for each read length with an offset, of each compartment
#'   (the rows of \code{data} of \code{results_choice}), and the columns \code{cutoff} (the offset),
#'   \code{frame_preference} (the percentage of the 5' ends in the 99 nt of the CDS windows that are in one frame,
#'   as the mean of the transcripts, for the frame with the highest mean), \code{read_length}, \code{gain_codons} (in
#'   the CDS windows, with the P-sites of this offset: the mean fraction of the positions in frame with P-sites,
#'   minus the fractions of the two other frames), \code{gain_new_codons} (the same, without the positions that the
#'   read lengths before it cover; see \code{\link{choose_readlengths}}), \code{all}, \code{max_inframe} and
#'   \code{max_coverage} (\code{TRUE} when the method of \code{readlength_choice_method} with this name selects the
#'   read length), \code{comp} (the compartment) and \code{pct_map} (the percentage of all reads of the BAM file that
#'   are of this compartment and read length, rounded to 4 digits). \code{""} when no read length has an offset.
#'   \code{offsets_df} and \code{rescue_all_rls} do not change it.}
#' }
#'
#' @section The smaller results file:
#' \code{_results_RiboseQC} has the list \code{res_all} with these changes: \code{read_stats$reads_pos1} and
#' \code{P_sites_stats} are \code{""}, and the profiles of \code{five_prime_bins}, \code{five_prime_subcodon},
#' \code{P_sites_bins} and \code{P_sites_subcodon} are summed over the transcripts: a list with one named numeric
#' vector (the column sums) for each compartment and read length. All other elements are the same.
#'
#' @section Files of the report:
#' With \code{create_report = TRUE}, after the last BAM file, \code{RiboseQC_analysis} makes the html report of the
#' \code{_results_RiboseQC} files with \code{\link{create_html_report}}. \code{<report>} is \code{report_file} (by
#' default, the first BAM file followed by \code{_RiboseQC_report.html}).
#' \describe{
#' \item{\code{<report>}}{The html report, with its plots in the file.}
#' \item{\code{<report>_plots/rds/}}{One file for each plot of the report, without a file extension, to read with
#'   \code{readRDS()}: a list of the plot (a ggplot, or a gtable for the heatmaps and the codon plots) and its width
#'   and height in inches. In the names, \code{<sample>} is the sample name (\code{sample_names}), \code{<comp>} the
#'   compartment, and the number is the section of the report:
#'   \itemize{
#'   \item \code{<sample>_1_readlocdist} and \code{all_samples_1_readlocdist}: the reads on each region (1).
#'   \item \code{<sample>_2_readlendist}: the reads of each read length (2).
#'   \item \code{<sample>_3_readloclendist1} and \code{<sample>_3_readloclendist2}: the reads by read length and
#'     region (3).
#'   \item \code{<sample>_<comp>_4_profiles_fivepr_metagene_<bins|subcodon>_<none|log2|zscore>}: the heatmaps of
#'     \code{five_prime_bins} and \code{five_prime_subcodon}, without a transformation, as \code{log2(x + 1)} or as
#'     z-scores (4.1). \code{<sample>_<comp>_4_profiles_fivepr_metagene_subcodon_<read length|all>}: the bar plots
#'     of \code{five_prime_subcodon}, for the read lengths (and \code{"all"}) with more than 100 5' ends.
#'   \item \code{<sample>_<comp>_4_selection5toP}: the frames of the 5' ends of each read length with an offset
#'     (4.2).
#'   \item \code{<sample>_<comp>_4_profiles_P_sites_metagene_...}: the same plots as in 4.1, of \code{P_sites_bins}
#'     and \code{P_sites_subcodon} (4.3); the bar plots for the read lengths with P-sites.
#'   \item \code{<sample>_<comp>_7_codonusage_positional_<type>_<read length|all>_<none|log2|zscore>}: the heatmaps
#'     of the codon tables, by position (7). \code{<type>} is \code{Codon_counts}, \code{P-sites_counts},
#'     \code{A-sites_counts} or \code{E-sites_counts} (\code{P_sites_percodon} and the others), or
#'     \code{P-sites_per_codon}, \code{A-sites_per_codon} or \code{E-sites_per_codon} (the \code{_ratio} tables).
#'   \item \code{<sample>_<comp>_8_codonusage_bulk_<type>_<read length|all>}: the codon tables summed over the
#'     positions (8).
#'   }
#'   The plots of sections 7 and 8 are only for \code{"all"}; with \code{extended_report = TRUE}, also for each read
#'   length of \code{Codon_counts}. Sections 4.2.2, 4.2.3, 5 and 6 are tables, without files.}
#' \item{\code{<report>_plots/pdf/}}{With \code{pdf_plots = TRUE}, one PDF file for each file of \code{rds/}, with its
#'   name followed by \code{.pdf} (see \code{\link{create_pdfs_from_rds_objects}}). The folder is always made.}
#' \item{\code{<report>_report_text_output.txt}}{The text output of rmarkdown while it makes the report: a log, not
#'   results.}
#' }
#'
#' @section The file of prepare_for_ORFquant:
#' \code{\link{prepare_for_ORFquant}} writes only \code{<dest_name>_for_ORFquant}, with the list \code{for_ORFquant}
#' as above. With the \code{_P_sites_offsets.tsv} of a run of \code{RiboseQC_analysis} and the same BAM file, it is the
#' same as the \code{_for_ORFquant} of that run, but without the reads on chromosomes that are not in the
#' annotation.
#'
#' @examples
#' # the smaller results file of the root sample of the example data
#' load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC"))
#' names(res_all)
#' res_all$read_stats$rld[, 1:10]
#' res_all$read_stats$positions
#' res_all$sequence_analysis[1:5]
#' # the root sample has no P-sites; the shoots sample has those of ChrC
#' load(system.file("ext_data", "shoots_results_RiboseQC", package = "RiboseQC"))
#' res_all$summary_P_sites
#' @seealso \code{\link{RiboseQC_analysis}}, \code{\link{create_html_report}}, \code{\link{prepare_for_ORFquant}},
#' \code{\link{prepare_annotation_files}}, \code{\link{calc_cutoffs_from_profiles}}, \code{\link{choose_readlengths}}
#' @name RiboseQC_output
#' @aliases res_all
NULL

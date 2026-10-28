# RiboseQC, a comprehensive Ribo-seq quality control tool
#
# Authors:
# Lorenzo Calviello (calviello.l.bio@gmail.com)
# Dominique Sydow (dominique.sydow@posteo.de)
# Dermott Harnett (Dermot.Harnett@mdc-berlin.de)
# Uwe Ohler (Uwe.Ohler@mdc-berlin.de)
#
# This software is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software. If not, see
# <http://www.gnu.org/licenses/>.

######################################################################################################
# Methods for RiboseQC report
######################################################################################################



#' Create the RiboseQC analysis report in html
#'
#' This function creates the RiboseQC html report based on the RiboseQC analysis files
#' generated with \code{RiboseQC_analysis}.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param input_files Character vector with full paths to data files
#' generated with \code{RiboseQC_analysis}.
#' Must be of same length as \code{input_sample_names}.
#'
#' @param input_sample_names Character vector containing input names (max. 5 characters per name).
#' Must be of same length as \code{input_files}.
#'
#' @param output_file String; full path to html report file.
#' @param extended creates a large html report including codon occupancy for each read length. Defaults to \code{FALSE}
#' @return The function saves the html report file with the file path \code{output_file} and
#' a folder containing all figures shown in the html report as RDS object files (located in the same directory as the html report).
#' It also saves the text output of the rendering in \code{<output_file>_report_text_output.txt}, and creates the folder
#' \code{<output_file>_plots/pdf/} for \code{create_pdfs_from_rds_objects}.
#'
#' @details This function creates the html report visualizing the RiboseQC analysis data. \cr \cr
#' Input are two character vectors of the same length: \cr \cr
#' a) \code{input_files}: vector of full paths to one or multiple input files
#' (RiboseQC analysis files generated with \code{RiboseQC_analysis}) and \cr \cr
#' b) \code{input_sample_names}: vector of corresponding names describing the file content in max. 5 characters
#' (these are used as names in the report). \cr \cr
#' For the report, a RMarkdown file is rendered as html document, saved as \code{output_file}. \cr \cr
#' Additionally, all figures in the report are saved as RDS object files
#' in an extra folder in the same directory as the report html file. \code{create_pdfs_from_rds_objects}
#' makes PDF figures of them (\code{RiboseQC_analysis} does this with \code{pdf_plots = TRUE}). \cr \cr
#' Example: \cr
#' \code{output_file <- "\\mydir\\myreport.html"} will generate
#' the html report \code{\\mydir\\myreport.html} and
#' the folder \code{\\mydir\\myreport.html_plots\\rds\\} for the RDS object files to be stored in.
#'
#' @seealso \code{\link{RiboseQC_output}} (the files of the report), \code{\link{plot_read_biotype_dist_1}},
#' \code{\link{plot_read_biotype_dist_2}}, \code{\link{plot_read_length_dist}}, \code{\link{plot_read_length_dist_by_biotype}},
#' \code{\link{plot_read_biotype_dist_by_length}},
#' \code{\link{get_metagene_data}},
#' \code{\link{plot_metagene_hm_rmd}}, \code{\link{plot_metagene_hm}},
#' \code{\link{plot_metagene_bar_rmd}}, \code{\link{plot_metagene_bar}},
#' \code{\link{plot_frame_dist_boxplot_rmd}}, \code{\link{plot_frame_dist_boxplot}},
#' \code{\link{get_rl_and_cutoffs}}, \code{\link{get_default_rl_selection}},
#' \code{\link{get_top50_mapping}}, \code{\link{get_top50_cds_genes}}, \code{\link{get_top50_all_genes}},
#' \code{\link{get_codon_usage_data}},
#' \code{\link{plot_codon_usage_positional_rmd}}, \code{\link{plot_codon_usage_positional}},
#' \code{\link{plot_codon_usage_bulk_rmd}}, \code{\link{plot_codon_usage_bulk}}
#' @examples
#' \dontrun{
#' sampnames = c("root","shoots")
#' create_html_report(input_files=paste(sampnames,"_results_RiboseQC",sep = ""), input_sample_names=sampnames,output_file =  'test_root_shoots.html',extended = FALSE)
#' }
#' @export

create_html_report <- function(input_files, input_sample_names, output_file,extended=FALSE){

    # get input and output file paths
    input_files <- paste(normalizePath(dirname(input_files)),basename(input_files),sep="/")
    output_file <- paste(normalizePath(dirname(output_file)),basename(output_file),sep="/")

    # get path to RMarkdown file (to be rendered)
    rmd_path <- paste(system.file(package="RiboseQC",mustWork = TRUE),"/rmd/riboseqc_template.Rmd",sep="")
    if(extended){
        rmd_path <- paste(system.file(package="RiboseQC",mustWork = TRUE),"/rmd/riboseqc_template_full.Rmd",sep="")
    }

    # get folder path for pdf figures
    output_fig_path <- paste(output_file,"_plots/", sep = "")

    # create folder for rds ojects and pdf figures
    dir.create(paste0(output_fig_path, "rds/"), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(output_fig_path, "pdf/"), recursive=TRUE, showWarnings=FALSE)
    sink(file = paste(output_file,"_report_text_output.txt",sep = ""))
    on.exit(sink())
    # render RMarkdown file > html report
    knitclean<-knitr::knit_meta(class=NULL, clean = TRUE)
    suppressWarnings(render(rmd_path,
           params = list(input_files = input_files,
                         input_sample_names = input_sample_names,
                         output_fig_path = output_fig_path),
           output_file = output_file,
           intermediates_dir = tempdir()))
    gici<-gc()
}


#' Generate PDF files from RDS object files
#'
#' This function generates figures as PDF files from RDS object files -
#' it's essentially a workaround for  difficulty integrating pdfs into
#' html reports.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param output_rds_path String; full path to output folder for RDS object files, ending with "/". Example: /my_path_to/rds/
#'
#' @return This function creates PDF files from RDS object files: for each RDS object file of a plot function in \code{output_rds_path},
#' the file \code{<name>.pdf} in the folder \code{pdf/} next to it (\code{<output_rds_path>../pdf/}). It returns \code{NULL}.
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC"))
#' output_rds_path <- file.path(tempdir(), "report_plots", "rds", "")
#' dir.create(output_rds_path, recursive = TRUE)
#' plot_read_length_dist(res_all$read_stats$rld, "root", output_rds_path)
#' create_pdfs_from_rds_objects(output_rds_path)
#' list.files(file.path(tempdir(), "report_plots", "pdf"))
#' @export

create_pdfs_from_rds_objects <- function(output_rds_path){
    rds_dir <- paste0(output_rds_path)
    rds_files <- list.files(rds_dir)

    for (rds_file in rds_files) {
        # read RDS data
        rds <- readRDS(paste0(output_rds_path, rds_file))

        # get plot data
        g <- rds[[1]]

        # get format info
        format <- rds[[2]]

        # make sure a pdf directory exists on same level as rds directory
        dir.create(paste0(output_rds_path, "../pdf/"), recursive=TRUE, showWarnings=FALSE)

        # some objects can be directly printed, others need as_ggplot (arrangeGrob-derived objects)
        if (class(g)[2] == "ggplot") {
            pdf(paste0(output_rds_path, "../pdf/", rds_file, ".pdf"),
                width=format[1], height=format[2])
            print(g)
            dev.off()
        }
        else {
            pdf(paste0(output_rds_path, "../pdf/", rds_file, ".pdf"),
                width=format[1], height=format[2])
            print(as_ggplot(g))
            dev.off()
        }
    }
}


#' Generate a list of R data objects
#'
#' This function generates a list of loaded RData files
#' to be used during RiboseQC report generation.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param input_files Named character vector of the paths of the 'results_RiboseQC' files generated by \code{RiboseQC_analysis};
#' the names are the sample names. \cr \cr
#' Example: \cr
#' \code{input_files <- c(sample1="//path//to//sample1", sample2="//path//to//sample2")}
#'
#' @return This function returns a list of loaded RData objects
#' that were generated by \code{RiboseQC_analysis}: the \code{res_all} of each file, with the names of \code{input_files}.
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' rdata_list <- generate_rdata_list(c(
#'   root = system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC")))
#' names(rdata_list)
#' @export

generate_rdata_list <- function(input_files){
    rdata_list <- list()
    for (i in seq_len(length(input_files))){
        load(input_files[i])
        rdata_list[[names(input_files)[i]]] <- res_all
    }
    return(rdata_list)
}


#' Plot read location distribution by biotype (and originating compartment)
#'
#' This function plots the read location distribution by biotype (and originating compartment)
#' for one input sample. \cr \cr
#' This plot is used in the RiboseQC report in section 1.1.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param pos \code{res_all$read_stats$positions} generated by \code{RiboseQC_analysis}. \cr \cr
#' \code{data.frame} containing the number of reads per biotype (rows) and originating compartment (columns).
#'
#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plot for the html report and returns \code{NULL}, invisibly.
#' With \code{output_rds_path}, it also saves the plot as the RDS object file \code{<output_rds_path><sample>_1_readlocdist}.
#'
#' @seealso \code{\link{create_html_report}}
#'
#' @examples
#' load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC"))
#' plot_read_biotype_dist_1(res_all$read_stats$positions, "root")
#' @export

plot_read_biotype_dist_1 <- function(pos, sample, output_rds_path=""){

    names(pos) <- gsub("reads_", "", names(pos))
    pos <- cbind(biotype=row.names(pos), pos)
    pos.m <- suppressMessages(melt(pos))
    pos.m$biotype <- factor(pos.m$biotype, levels=row.names(pos))
    pos.m$variable <- factor(pos.m$variable, levels=names(pos))
    head(pos.m)

    g2 <- ggplot(data = pos.m, aes(x=biotype, y=value, fill=variable)) +
        geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
        labs(x="", y="read count\n") +
        scale_fill_brewer(palette = "Set2") +
        theme_minimal() +
        theme(plot.margin = margin(.5, .5, .5, .5, "cm"),
              panel.background=element_blank(),
              axis.ticks=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=14),
              axis.text.y=element_text(size=14),
              axis.title=element_text(size=14),
              legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title=element_blank())

    # print to HTML file
    print(g2)

    cat("\n\n")  # necessary for RMarkdown rendering

    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_1_readlocdist")
        saveRDS(list(g2, c(10,5)), file = file_path,
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}


#' Plot read location distribution by originating compartment (and biotype)
#'
#' This function plots the read location distribution by originating compartment (and biotype)
#' for all input samples. \cr \cr
#' This plot is used in the RiboseQC report in section 1.2.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plot for the html report and returns \code{NULL}, invisibly.
#' With \code{output_rds_path}, it also saves the plot as the RDS object file \code{<output_rds_path>all_samples_1_readlocdist}.
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC"))
#' plot_read_biotype_dist_2(list(root = res_all))
#' @export

plot_read_biotype_dist_2 <- function(rdata_list, output_rds_path=""){
    levels_bio <- NULL
    levels_comp <- NULL
    levels_sample <- names(rdata_list)

    data <- list()
    for (i in names(rdata_list)){
        res_all <- rdata_list[[i]]
        rs <- res_all$read_stats$reads_summary
        levels_bio <- row.names(rs[[1]])
        levels_comp <- names(rs)
        data[[i]] <- do.call(cbind, lapply(rs, function(x) rowSums(as.data.frame(x))))
    }

    data.m <- suppressMessages(melt(data))
    data.m$Var1 <- factor(data.m$Var1, levels=levels_bio)
    data.m$Var2 <- factor(data.m$Var2, levels=levels_comp)
    data.m$L1 <- factor(data.m$L1, levels=levels_sample)

    # read count
    g1 <- ggplot(data=data.m, aes(x=L1, y=value, fill=Var1)) +
        geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
        scale_fill_viridis(discrete=TRUE) +
        labs(x="\nsample", y="read count\n") +
        facet_grid(Var2~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              panel.background=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=14),
              axis.text.y=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              strip.text.y = element_text(size=14)) +
        guides(fill=guide_legend(nrow=1))

    # read count percentage
    g2 <- ggplot(data=data.m, aes(x=L1, y=value, fill=Var1)) +
        geom_bar(stat="identity", position=position_fill(reverse=TRUE)) +
        scale_fill_viridis(discrete=TRUE) +
        labs(x="\nsample", y="read count fraction\n") +
        facet_grid(Var2~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              panel.background=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=14),
              axis.text.y=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              strip.text.y = element_text(size=14)) +
        guides(fill=guide_legend(nrow=1))

    gg <- ggarrange(g1, g2, ncol = 2, common.legend=TRUE, legend="bottom") +
        theme(plot.margin = margin(.5, .5, .5, .5, "cm"))

    # print to HTML file
    print(gg)

    cat("\n\n")  # necessary for RMarkdown rendering

    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, "all_samples_1_readlocdist")
        saveRDS(list(gg, c(10,7)), file = file_path,
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}



#' Plot read length distribution
#'
#' This function plots the read length distribution per originating compartment
#' for one input sample. \cr \cr
#' This plot is used in the RiboseQC report in section 2.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rld \code{res_all$read_stats$rld} generated by \code{RiboseQC_analysis} \cr \cr
#' \code{data.frame} containing the number of reads per originating compartment (rows) and read lengths (columns).
#'
#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plot for the html report and returns \code{NULL}, invisibly.
#' With \code{output_rds_path}, it also saves the plot as the RDS object file \code{<output_rds_path><sample>_2_readlendist}.
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC"))
#' plot_read_length_dist(res_all$read_stats$rld, "root")
#' @export

plot_read_length_dist <- function(rld, sample, output_rds_path=""){
    lw <- 1.5 # line width

    rld <- as.matrix(rld)
    colnames(rld) <- as.integer(gsub("reads_", "", colnames(rld)))

    levels_comp <- rownames(rld)

    # read count
    rld.m <- suppressMessages(melt(rld))
    rld.m$Var1 <- factor(rld.m$Var1, levels=levels_comp)
    rld.m$Var2 <- as.numeric(rld.m$Var2)
    colssi<-alpha(c("blue","dark red","forestgreen",rainbow(length(unique(rld.m$Var1))-3)),alpha = .7)



    g1 <- ggplot(data=rld.m, aes(x=Var2, y=value)) +
        geom_line(aes(colour=Var1, group=Var1), linewidth=lw) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        labs(x="\nread length", y="read count\n") +
        scale_color_manual(values = colssi) +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank())

    # read count fraction
    rld_norm <- rld/rowSums(rld)
    rld_norm.m <- suppressMessages(melt(rld_norm))
    rld_norm.m$Var1 <- factor(rld_norm.m$Var1, levels=levels_comp)
    colssi<-alpha(c("blue","dark red","forestgreen",rainbow(length(unique(rld.m$Var1))-3)),alpha = .7)

    g2 <- ggplot(data=rld_norm.m, aes(x=Var2, y=value)) +
        geom_line(aes(colour=Var1, group=Var1), linewidth=lw) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        labs(x="\nread length", y="read count fraction\n") +
        scale_color_manual(values = colssi) +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank())

    gg <- ggarrange(g1, g2, ncol = 2, common.legend=TRUE, legend="bottom") +
        theme(plot.margin = margin(.5, .5, .5, .5, "cm"))

    # print to HTML file
    print(gg)

    cat("\n\n")  # necessary for RMarkdown rendering

    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_2_readlendist")
        saveRDS(list(gg, c(10,5)), file = file_path,
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }

}



#' Plot read length and location distribution
#' (distribution per biotype)
#'
#' This funtion plots the read length distribution for each originating compartment as distribtion per biotype
#' for one input sample (displayed as read count and as read count fraction). \cr \cr
#' This plot is used in the RiboseQC report in section 3.1.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param reads_summary \code{res_all$read_stats$reads_summary} generated by \code{RiboseQC_analysis} \cr \cr
#' List of DataFrames:
#' one DataFrame for each originating compartment,
#' each DataFrame contains read counts per biotype (rows) and read lengths (columns).
#'
#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plot for the html report and returns \code{NULL}, invisibly.
#' With \code{output_rds_path}, it also saves the plot as the RDS object file \code{<output_rds_path><sample>_3_readloclendist1}.
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC"))
#' plot_read_length_dist_by_biotype(res_all$read_stats$reads_summary, "root")
#' @export

plot_read_length_dist_by_biotype <- function(reads_summary, sample, output_rds_path=""){
    lw <- 1.5 # line width

    levels_bio <- row.names(reads_summary[[1]])
    levels_comp <- names(reads_summary)

    rld <- lapply(reads_summary, function(x) as.matrix(x))

    # read count
    rld.m <- suppressMessages(melt(rld))
    rld.m$Var1 <- factor(rld.m$Var1, levels=levels_bio)
    rld.m$Var2 <- as.integer(gsub("reads_", "", rld.m$Var2))
    rld.m$L1 <- factor(rld.m$L1, levels=levels_comp)

    g1 <- ggplot(data=rld.m, aes(x=Var2, y=value)) +
        geom_line(aes(colour=Var1, group=Var1), linewidth=lw) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        labs(x="\nread length", y="read count\n") +
        scale_color_viridis(discrete=TRUE,alpha = .8) +
        facet_grid(L1~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              strip.text.y = element_text(size=14)) +
        guides(colour=guide_legend(nrow=1))

    # read count fraction
    rld_norm <- list()
    for (i in names(rld)){
        rld_norm[[i]] <- rld[[i]]/(rowSums(rld[[i]]))
        rld_norm[[i]][is.na(rld_norm[[i]])] = 0 # replace NaN fields with 0
    }
    rld_norm.m <- suppressMessages(melt(rld_norm))
    rld_norm.m$Var1 <- factor(rld_norm.m$Var1, levels=levels_bio)
    rld_norm.m$Var2 <- as.integer(gsub("reads_", "", rld_norm.m$Var2))
    rld_norm.m$L1 <- factor(rld_norm.m$L1, levels=levels_comp)

    g2 <- ggplot(data=rld_norm.m, aes(x=Var2, y=value)) +
        geom_line(aes(colour=Var1, group=Var1), linewidth=lw) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        labs(x="\nread length", y="read count fraction\n") +
        scale_color_viridis(discrete=TRUE,alpha = .8) +
        facet_grid(L1~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              strip.text.y = element_text(size=14)) +
        guides(colour=guide_legend(nrow=1))

    gg <- ggarrange(g1, g2, ncol = 2, common.legend=TRUE, legend="bottom") +
        theme(plot.margin = margin(.5, .5, .5, .5, "cm"))

    # print to HTML file
    plot(gg)

    cat("\n\n")  # necessary for RMarkdown rendering

    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_3_readloclendist1")
        saveRDS(list(gg, c(10,8)), file = file_path,
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}



#' Plot read length and location distribution
#' (distribution per read length)
#'
#' This funtion plots the biotype distribution for each originating compartment as distribtion per read length
#' for one input sample (displayed as read count and as read count fraction). \cr \cr
#' This plot is used in the RiboseQC report in section 3.2.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param reads_summary \code{res_all$read_stats$reads_summary} generated by \code{RiboseQC_analysis} \cr \cr
#' List of DataFrames:
#' one DataFrame for each originating compartment,
#' each DataFrame contains read counts per biotype (rows) and read lengths (columns).

#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plot for the html report and returns \code{NULL}, invisibly.
#' With \code{output_rds_path}, it also saves the plot as the RDS object file \code{<output_rds_path><sample>_3_readloclendist2}.
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC"))
#' plot_read_biotype_dist_by_length(res_all$read_stats$reads_summary, "root")
#' @export

plot_read_biotype_dist_by_length <- function(reads_summary, sample, output_rds_path=""){
    levels_bio <- row.names(reads_summary[[1]])
    levels_comp <- names(reads_summary)

    rld <- lapply(reads_summary, function(x) as.matrix(x))

    # read count
    rld.m <- suppressMessages(melt(rld))
    rld.m$Var1 <- factor(rld.m$Var1, levels=levels_bio)
    rld.m$Var2 <- as.integer(gsub("reads_", "", rld.m$Var2))
    rld.m$L1 <- factor(rld.m$L1, levels=levels_comp)


    # read count
    g1 <- ggplot(data=rld.m, aes(x=Var2, y=value, fill=Var1)) +
        geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        scale_fill_viridis(discrete=TRUE) +
        labs(x="\nread length (nt)", y="read count\n") +
        facet_grid(L1~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              panel.background=element_blank(),
              axis.text.x=element_text(angle=90, vjust=0.5, size=14),
              axis.text.y=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              strip.text.y = element_text(size=14)) +
        guides(fill=guide_legend(nrow=1))

    # read count percentage
    g2 <- ggplot(data=rld.m, aes(x=Var2, y=value, fill=Var1)) +
        geom_bar(stat="identity", position=position_fill(reverse=TRUE)) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        scale_fill_viridis(discrete=TRUE) +
        labs(x="\nread length (nt)", y="read count fraction\n") +
        facet_grid(L1~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              panel.background=element_blank(),
              axis.text.x=element_text(angle = 90, vjust=0.5, size=14),
              axis.text.y=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              legend.position = "bottom",
              strip.text.y = element_text(size=14)) +
        guides(fill=guide_legend(nrow=1))

    gg <- ggarrange(g1, g2, ncol = 2, common.legend=TRUE, legend="bottom") +
        theme(plot.margin=margin(.5, .5, .5, .5, "cm"))

    # print to HTML file
    print(gg)

    cat("\n\n")  # necessary for RMarkdown rendering

    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_3_readloclendist2")
        saveRDS(list(gg, c(10,8)), file = file_path,
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}


#' Get 5'/P-site profile data for metagene analysis
#'
#' This function processes profile data generated by \code{RiboseQC_analysis}. \cr \cr
#' This data is used as input in \code{plot_metagene_hm} to generate plot for the RiboseQC report in section 4.1/4.3.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param data Object (list of lists) generated by \code{RiboseQC_analysis}: \code{res_all} of the 'results_RiboseQC' file.
#'
#' @param profile_type String; \code{"profiles_fivepr"} or \code{"profiles_P_sites"}:
#' the 5' or P-site profile data generated by \code{RiboseQC_analysis}: \cr \cr
#' \code{res_all$profiles_fivepr} or \cr
#' \code{res_all$profiles_P_sites} \cr \cr
#' Consists of numeric vectors (the profiles summed over the transcripts; the DataFrames of the 'results_RiboseQC_all' file
#' do not work here) each
#' containing counts of 5' or P-site profiles,
#' calculated for different resolution types (see parameter \code{res}),
#' originating compartments (see parameter \code{comp}), and
#' read lengths per input sample. \cr \cr
#' Example to access a profile: \cr
#' \code{res_all$profiles_fivepr[[res]][[comp]][[read_length]]} or \cr
#' \code{res_all$profiles_P_sites[[res]][[comp]][[read_length]]}
#'
#' @param comp String for originating compartment \cr \cr
#' Check for available originating compartments in the data set using:
#' \code{names(res_all$profiles_fivepr$five_prime_subcodon)}
#'
#' @param res Resolution (subcodon or bins)  \cr \cr
#' \bold{Subcodon} resolution: \cr
#' \code{five_prime_subcodon} \cr
#' (in order to call \code{res_all$profiles_fivepr$five_prime_subcodon}) or \cr
#' \code{P_sites_subcodon} \cr
#' (in order to call \code{res_all$profiles_P_sites$P_sites_subcodon}) \cr \cr
#' Read coverage for
#' the first 25nt after the transcription start site (TSS),
#' 25nt before the start codon and the first 33nt of the CDS (from the start codon),
#' 33nt from the middle of the CDS,
#' the last 33nt of the CDS and 25nt after the stop codon,
#' and the last 25nt before the transcription end site (TES). \cr \cr
#' \bold{Bins}: \cr
#' \code{five_prime_bins} \cr
#' (in order to call \code{res_all$profiles_fivepr$five_prime_bins}) or \cr
#' \code{P_sites_bins} \cr
#' (in order to call \code{res_all$profiles_P_sites$P_sites_bins}) \cr \cr
#' Read coverage for
#' 50 bins between TSS and start codon, 100 bins for the CDS, and 50 after stop codon to TES.
#'
#' @return This function returns a list with the profile data:
#' \itemize{
#' \item \code{data_bar}: for each read length and for all read lengths summarized (\code{"all"}), a data.frame with
#' the profile (\code{value}, \code{pos} and \code{frame}), for \code{plot_metagene_bar},
#' \item \code{data_hm_single}: the profiles of all read lengths individually, for the heatmaps, and
#' \item \code{data_hm_all}: the profile of all read lengths summarized (\code{"all"}), for the heatmaps,
#' }
#' with different scaling: no scaling (\code{none}), log2 scaling (\code{log2}), and z scoring (\code{zscore}),
#' accessable via e.g. \code{profile_data$data_hm_single$none} or \code{profile_data$data_hm_all$zscore}. \cr \cr
#' \code{profile_type}, \code{res} and \code{comp} save the arguments for later use during plotting.
#'
#' @seealso \code{\link{create_html_report}}
#' @export

get_metagene_data <- function(data, profile_type, res, comp){

    # test input parameters for existence
    if((profile_type != "profiles_fivepr") & (profile_type != "profiles_P_sites")){
        stop("Variable profile_type takes only: \"profiles_fivepr\" or \"profiles_P_sites\"")
    }
    if(profile_type == "profiles_fivepr") {
        if((res != "five_prime_bins") & (res != "five_prime_subcodon")){
            stop("Variable res takes only: \"five_prime_bins\" or \"five_prime_subcodon\"")
        }
    }
    if(profile_type == "profiles_P_sites") {
        if((res != "P_sites_bins") & (res != "P_sites_subcodon")){
            stop("Variable res takes only: \"P_sites_bins\" or \"P_sites_subcodon\"")
        }
    }
    if(!(comp %in% names(data[[profile_type]][[res]]))) {
        stop("Variable comp takes only: ", paste(names(data[[profile_type]][[res]]), collapse=", "))
    }

    # get metagene data
    signal <- data[[profile_type]][[res]][[comp]]
    # get available read lengths
    rl_ok <- names(signal)
    # sort read lengths in ascending order and "all" at the beginning
    rl_ok <- c(rl_ok[1], as.character(sort(as.integer(rl_ok[-1], decreasing=FALSE))))

    # prepare BARPLOTS data

    data_bar = list()

    for (rl in rl_ok) {
        # get read length data and format
        x <- signal[[rl]]
        x.m <- suppressMessages(melt(x))
        x.m$pos <- seq_len(length(x))
        x.m$frame <- paste0("frame ", rep(c(3,1,2),67)[seq_along(x)])

        # add read length data to list
        data_bar[[rl]] <- x.m
    }

    # prepare HEATMAP data

    # we want to make two separate heatmaps for
    # a) single read lenghts (single) and
    # b) summarized read lengths (all).

    # single read lengths
    rl_ok_single <- rl_ok[-1]
    rl_ok_single <- c(as.character(sort(as.integer(rl_ok_single), decreasing = FALSE)))
    data_none <- do.call(signal[rl_ok_single], what=rbind)
    rownames(data_none) <- rl_ok_single
    colnames(data_none) <- seq_len(dim(data_none)[2])
    data_log2 <- log2(data_none+1)
    data_zscore <- scale(data_none, center=TRUE, scale=TRUE) # scale centers (true) and scales (true) the columns of a numeric matrix, we want to scale per read length > transpose data_none to have rl as column
    data <- list(none=suppressMessages(melt(t(data_none))), log2=suppressMessages(melt(t(data_log2))), zscore=suppressMessages(melt(t(data_zscore))))
    data$none$Var2 <- factor(data$none$Var2,levels=rl_ok_single)
    data$log2$Var2 <- factor(data$log2$Var2,levels=rl_ok_single)
    data$zscore$Var2 <- factor(data$zscore$Var2,levels=rl_ok_single)
    data_hm_single <- data

    # read length "all"
    rl_ok_all <- rl_ok[1]
    data_none <- do.call(signal[rl_ok_all], what=rbind)
    rownames(data_none) <- names(signal[rl_ok_all])
    colnames(data_none) <- seq_len(dim(data_none)[2])
    data_log2 <- log2(data_none+1)
    #it needs the transpose as it's a 1-row matrix
    data_zscore <- t(scale(t(data_none), center=TRUE, scale=TRUE))
    data <- list(none=suppressMessages(melt(t(data_none))), log2=suppressMessages(melt(t(data_log2))), zscore=suppressMessages(melt(t(data_zscore))))
    data$none$Var2 <- data$none$Var2
    data$log2$Var2 <- data$log2$Var2
    data$zscore$Var2 <- data$zscore$Var2
    data_hm_all <- data

    # save barplot, heatmap and meta data!
    metagene_data <- list(data_bar=data_bar,
                          data_hm_single=data_hm_single,
                          data_hm_all=data_hm_all,
                          profile_type=profile_type,
                          res=res,
                          comp=comp)

    return(metagene_data)
}



#' Plot 5'/P-site profiles as heatmaps within the RMarkdown document
#'
#' This function generates iteratively all heatmap plots
#' (iteration over originating compartments, resolution types, and scaling methods),
#' and for subcodon resolution without scaling the barplots of \code{plot_metagene_bar_rmd}. \cr \cr
#' These plots are displayed in the RiboseQC report in section 4.1/4.3.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param data Object (list of lists) generated by \code{RiboseQC_analysis}: \code{res_all} of the 'results_RiboseQC' file.
#'
#' @param profile_type String; \code{"profiles_fivepr"} or \code{"profiles_P_sites"}:
#' the 5' or P-site profile data generated by \code{RiboseQC_analysis}: \cr \cr
#' \code{res_all$profiles_fivepr} or \cr
#' \code{res_all$profiles_P_sites} \cr \cr
#' Consists of numeric vectors (the profiles summed over the transcripts; the DataFrames of the 'results_RiboseQC_all' file
#' do not work here) each
#' containing counts of 5' or P-site profiles,
#' calculated for different resolution types (see parameter \code{res}),
#' originating compartments (see parameter \code{comp}), and
#' read lengths per input sample. \cr \cr
#' Example to access a profile: \cr
#' \code{res_all$profiles_fivepr[[res]][[comp]][[read_length]]} or \cr
#' \code{res_all$profiles_P_sites[[res]][[comp]][[read_length]]}
#'
#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints iteratively all 5' or P-site profile plots
#' for the html report and returns \code{NULL}. With \code{output_rds_path}, it also saves the same plots as RDS object files.
#'
#' @seealso \code{\link{create_html_report}}, \code{\link{get_metagene_data}}, \code{\link{plot_metagene_hm}}
#' @examples
#' load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC"))
#' plot_metagene_hm_rmd(res_all, "profiles_fivepr", "root")
#' @export

plot_metagene_hm_rmd <- function(data, profile_type, sample="", output_rds_path=""){

    # split data by resolution (subcodon vs. bins)
    profile <- data[[profile_type]]
    ress <- rev(names(profile))
    ress <- c(ress[grepl("subcodon", ress)], ress[grepl("bins", ress)])

    comps <- names(profile[[1]])

    for (comp in comps){
        # cat commands necessary for RMarkdown rendering
        cat("\n\n")
        cat("#### ", comp, " {.tabset} \n\n")
        cat("Select a *resolution* (subcodon or bins):\n")
        cat("In case of *subcodon* resolution, read coverage is shown for
            the first 25nt after the transcription start site (TSS),
            25nt before the start codon and the first 33nt of the CDS (from the start codon),
            33nt from the middle of the CDS,
            the last 33nt of the CDS and 25nt after the stop codon,
            and the last 25nt before the transcription end site (TES).\n")
        cat("In case of *bins*, read coverage is shown for 50 bins between TSS and start codon (5'UTR),
            100 bins for the CDS, and 50 after stop codon (3'UTR). \n\n")
        cat("Select a *scaling* method: *none* (no scaling), *log2* scaling, and *z-score* scaling.\n\n")
        cat("**Note**: Select *subcodon / none* or *bins / none* scaling to display the read coverage for each read length separately.\n\n")

        for(res in ress){

            # get heatmap data
            metagene_data <- get_metagene_data(data, profile_type, res, comp)

            # get different scaling factors (none, log2, z-score)
            scals <- names(metagene_data$data_hm_single)

            for (scal in scals){
                cat("\n\n")  # necessary for RMarkdown rendering

                # get resolution-scale combo, e.g. subcodon / log2
                cat("##### ", strsplit(res, split="_")[[1]][3], "/", scal, " {.tabset} \n\n")  # necessary for RMarkdown rendering

                # plot heatmap (displaying single read lengths and summarized read lengths (all))
                plot_metagene_hm(metagene_data, scal, sample, output_rds_path)
            }
        }
    }
}


#' Plot 5'/P-site profiles (for all read lengths) as heatmap
#'
#' This function plots a 5' or P-site profiles for all read lengths as heatmap
#' (for a specific originating compartment, resolution type, and scaling method).
#' With \code{scal = "none"} and subcodon resolution, it also plots the barplots of \code{plot_metagene_bar_rmd}. \cr \cr
#' This plot is used in the RiboseQC report in section 4.1/4.3.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param metagene_data Profile data for a specific originating compartment and resolution type,
#' generated using \code{\link{get_metagene_data}}.
#'
#' @param scal Scaling method:
#' no scaling (\code{none}), log2 scaling (\code{log2}), or z scoring (\code{zscore}).
#'
#' @param sample String; sample name (selected from the input names
#' iven in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plot for the html report and returns \code{NULL}, invisibly.
#' With \code{output_rds_path}, it also saves the plot as the RDS object file
#' \code{<output_rds_path><sample>_<comp>_4_<profile_type>_metagene_<subcodon or bins>_<scal>}.
#'
#' @seealso \code{\link{create_html_report}}
#' @export

plot_metagene_hm <- function(metagene_data, scal, sample="", output_rds_path="") {

    # heatmap data
    data_hm_single <- metagene_data$data_hm_single
    data_hm_all <- metagene_data$data_hm_all

    # metadata
    profile_type <- metagene_data$profile_type
    res <- metagene_data$res
    comp <- metagene_data$comp

    # set subcodon/binned resolution-specific parameters
    if (grepl("subcodon", res)) {
        res_breaks <- c(1, 26, 51, 84, 117, 150, 176, 200)
        res_labels <- c("TSS", "", "start\ncodon", "", "", "stop\ncodon", "", "TES")
        res_x_lab <- "\nposition (nucleotide resolution)"

    } else {
        res_breaks <- c(1, 51, 150, 200)
        res_labels <- c("TSS", "start\ncodon", "stop\ncodon", "TES")
        res_x_lab <- "\nbins"
    }

    # heatmap for summarized read lengths (all)
    g1 <- ggplot(data_hm_all[[scal]], aes(x=Var1, y=Var2)) +
        geom_tile(aes(fill=value)) +
        scale_fill_gradient(low="white", high="steelblue", na.value="gray80") +
        scale_x_continuous(breaks=res_breaks, labels=res_labels) +
        scale_y_discrete(expand=c(0, 0)) +
        theme(axis.text=element_text(size=18),
              axis.title.y=element_text(size=18, colour="white"),
              axis.ticks.y=element_blank(),
              panel.background=element_blank(),
              axis.title.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.title=element_text(size=18),
              legend.text=element_text(size=14)) +
        labs(colour="r", y="lg\n", x=res_x_lab) +
        geom_vline(xintercept=c(51,150), lty=3, color="black") +
        geom_point(data=data_hm_all[[scal]], aes(size="NA"), shape=NA, colour="gray80") +
        guides(fill=guide_legend(title=paste0("read count (scale: ", scal, ")"), nrow=1),
               size=guide_legend(title=""))

    data_hm_single[[scal]]$Var2<-as.numeric(as.character(data_hm_single[[scal]]$Var2))
    # heatmap for single read lenghts (single)
    g2 <- ggplot(data_hm_single[[scal]], aes(x=Var1, y=Var2)) +
        geom_tile(aes(fill=value)) +
        scale_fill_gradient(low="white", high="steelblue", na.value="gray80") +
        scale_x_continuous(breaks=res_breaks, labels=res_labels) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=18),
              panel.background=element_blank(),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.title=element_text(size=18),
              legend.text=element_text(size=14,angle = 45)) +
        labs(colour="r", y="read length (nt)\n", x=res_x_lab) +
        geom_vline(xintercept=c(51,150), lty=3, color="black") +
        geom_point(data=data_hm_single[[scal]], aes(size="NA"), shape=NA, colour="gray80") +
        labs(fill=paste0("read count (scale: ", scal, ")")) + guides(size=guide_legend(title=""))

    # zscore: the plot with other colours below
    if(scal!="zscore"){gg <- arrangeGrob(g1, g2, ncol=1, heights = c(1, 6))}

    suppressMessages(if(scal=="zscore"){gg<-arrangeGrob(g1 + scale_fill_gradient2(low = "blue",mid = "white",high = "firebrick",na.value = "gray80",guide = "colourbar"),
                    g2 + scale_fill_gradient2(low = "blue",mid = "white",high = "firebrick",na.value = "gray80",guide = "colourbar"),
                    ncol=1, heights = c(1, 6))})



    cat("\n\n")  # necessary for RMarkdown rendering

    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path,
                           sample, "_",
                           comp, "_",
                           "4_", profile_type, "_",
                           "metagene_", unlist(strsplit(res, "_"))[3], "_",
                           scal)
        saveRDS(list(gg, c(12,10)), file = file_path,
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }

    # print to HTML file
    print(as_ggplot(gg))

    # generate per-read-length metagene plots if scaling is none
    if(scal=="none" & grepl("subcodon", res)){
        plot_metagene_bar_rmd(metagene_data, sample, output_rds_path)
    }
}


#' Plot 5'/P-site profile per read length as barplot in R Markdown
#'
#' This function plots the 5' or P-site profiles as barplots with \code{plot_metagene_bar}
#' (for a specific originating compartment and resolution type): one for each read length and for all read lengths
#' summarized (\code{"all"}) with a profile, with a sum of more than 100 for 5' profiles and more than 0 for P-site profiles. \cr \cr
#' This plot is used in the RiboseQC report in section 4.1/4.3.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param metagene_data Profile data for a specific originating compartment and resolution type,
#' generated using \code{\link{get_metagene_data}}.
#'
#'
#' @param sample String; sample name (selected from the input names given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function. Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plots for the html report and returns \code{NULL}.
#' With \code{output_rds_path}, it also saves the plots as RDS object files (see \code{plot_metagene_bar}).
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' load(system.file("ext_data", "root_results_RiboseQC", package = "RiboseQC"))
#' metagene_data <- get_metagene_data(res_all, "profiles_fivepr", "five_prime_subcodon", "nucl")
#' plot_metagene_bar_rmd(metagene_data, "root")
#' @export

plot_metagene_bar_rmd <- function(metagene_data, sample="", output_rds_path="") {

    data_bar <- metagene_data$data_bar
    profile_type <- metagene_data$profile_type
    rls <- names(data_bar) # read lengths

    # first, check which read lengths you really want to display:
    # only keep read lengths with a cutoff coverage
    rls_selected = rls

    # do this only for the 5' profiles
    # (for P-site profiles, we want to see everything that is there, except empty read lengths)
    # profile_type 1 refers to the 5' profiles!
    if(profile_type=="profiles_fivepr"){
        # set cutoff (coverage per read length)
        cutoff <- 100
        # calculate coverage per read length
        rls_cov <- unlist(lapply(lapply(data_bar, "[", , 1), sum))[rls]
        # select read lengths according to cutoff
        rls_selected <- names(rls_cov[rls_cov > cutoff])
    }
    else{
        # set cutoff (coverage per read length)
        cutoff <- 0
        # calculate coverage per read length
        rls_cov <- unlist(lapply(lapply(data_bar, "[", , 1), sum))[rls]
        # select read lengths according to cutoff
        rls_selected <- names(rls_cov[rls_cov > cutoff])
    }

    for (rl in rls_selected){
        cat("\n\n")
        cat("###### ", rl, " \n\n")

        plot_metagene_bar(metagene_data, rl, sample, output_rds_path)
    }
}


#' Plot 5'/P-site profile per read length as barplot
#'
#' This function plots a 5' or P-site profile as barplot
#' (for a specific originating compartment, resolution type, and read length). \cr \cr
#' This plot is used in the RiboseQC report in section 4.1/4.3.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param metagene_data Profile data for a specific originating compartment and resolution type,
#' generated using \code{\link{get_metagene_data}}.
#'
#' @param rl String for read length: \code{"all"} (all read lengths) or a read length, for example \code{"30"}. \cr
#' Check for available read lengths using: \code{names(metagene_data$data_bar)}
#'
#'
#' @param sample String; sample name (selected from the input names
#' iven in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plot for the html report and returns \code{NULL}, invisibly.
#' With \code{output_rds_path}, it also saves the plot as the RDS object file
#' \code{<output_rds_path><sample>_<comp>_4_<profile_type>_metagene_<subcodon or bins>_<rl>}.
#'
#' @seealso \code{\link{create_html_report}}
#' @export

plot_metagene_bar <- function(metagene_data, rl, sample="", output_rds_path="") {

    # bar data
    data_bar <- metagene_data$data_bar
    # applied to selected read length
    x.m <- data_bar[[rl]]

    # metadata
    profile_type <- metagene_data$profile_type
    res <- metagene_data$res
    comp <- metagene_data$comp

    # set subcodon/binned resolution-specific parameters
    if(grepl("subcodon", res)) {
        res_breaks <- c(1, 26, 51, 84, 117, 150, 176, 200)
        res_labels <- c("TSS", "", "start\ncodon", "", "", "stop\ncodon", "", "TES")
        res_x_lab <- "\nposition (nucleotide resolution)"

    } else {
        res_breaks <- c(1, 51, 150, 200)
        res_labels <- c("TSS", "start\ncodon", "stop\ncodon", "TES")
        res_x_lab <- "\nbins"
    }

    # we need to treat subcodon and binned data differently:
    # subcodon data has frame information, therefore frames should be colored accordingly
    # binned data has NO frame information, therefore colour everything in black, please, dear computer

    # plot data
    if(grepl("subcodon", res)){

        x.m$frame<-as.character(x.m$frame)
        x.m$frame[x.m$frame=="frame 2"]<-"Frame 0"
        x.m$frame[x.m$frame=="frame 1"]<-"Frame 2"
        x.m$frame[x.m$frame=="frame 3"]<-"Frame 1"
        x.m$frame<-gsub(x.m$frame,pattern = "F",replacement = "f")
        x.m$frame<-factor(x.m$frame,levels=c("frame 0","frame 1","frame 2"))

        g3 <- ggplot(data=x.m, aes(x=pos, xend=pos, y=0, yend=value)) +
            geom_segment(aes(colour=frame, group=frame)) +
            scale_color_manual(values = c("red","forestgreen","blue")) +
            scale_x_continuous(breaks=res_breaks, labels=res_labels) +
            labs(x=res_x_lab, y="read count\n") +
            theme_minimal() +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 18),
                  legend.title = element_blank(),
                  axis.text = element_text(size = 18),
                  axis.title = element_text(size = 18))
    } else {
        g3 <- ggplot(data=x.m, aes(x=pos, xend=pos, y=0, yend=value)) +
            geom_segment(aes(group=frame), colour="black") +
            scale_x_continuous(breaks=res_breaks, labels=res_labels) +
            labs(x=res_x_lab, y="read count\n") +
            theme_minimal() +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 18),
                  legend.title = element_blank(),
                  axis.text = element_text(size = 18),
                  axis.title = element_text(size = 18))
    }

    # print to HTML file
    print(g3)

    cat("\n\n")  # necessary for RMarkdown rendering

    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_",
                           comp, "_",
                           "4_", profile_type, "_",
                           "metagene_",
                           unlist(strsplit(res, "_"))[3], "_",
                           rl)
        saveRDS(list(g3, c(12,10)), file = file_path,
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}


#' Plot frame coverage of 5'-site profiles within the RMarkdown document
#'
#' This function generates iteratively all frame coverage boxplots
#' (iteration over originating compartments). \cr \cr
#' These plots are displayed in the RiboseQC report in section 4.2.1.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param analysis_frame_cutoff Object containing statistics on frame coverage
#' (per originating compartment and read length) generated by \code{RiboseQC_analysis}. \cr \cr
#' Example: \cr
#' \code{res_all$selection_cutoffs$analysis_frame_cutoff}
#'
#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function integrates plots in the html report, for each compartment with data, and returns \code{NULL}.
#' With \code{output_rds_path}, it also saves the plots as RDS object files (see \code{plot_frame_dist_boxplot}).
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' # the shoots sample: no read length of the root sample has an offset
#' load(system.file("ext_data", "shoots_results_RiboseQC", package = "RiboseQC"))
#' plot_frame_dist_boxplot_rmd(res_all$selection_cutoffs$analysis_frame_cutoff, "shoots")
#' @export

plot_frame_dist_boxplot_rmd <- function(analysis_frame_cutoff, sample="", output_rds_path=""){
    comps <- names(analysis_frame_cutoff)
    # only compartments with data
    comps <- comps[unlist(lapply(comps, function(x) length(analysis_frame_cutoff[[x]])))!=0]

    for (comp in comps){
        cat("\n\n")  # necessary for RMarkdown rendering
        cat("##### ", comp, " {.tabset} \n\n")  # necessary for RMarkdown rendering

        plot_frame_dist_boxplot(analysis_frame_cutoff, comp, sample, output_rds_path)
    }
}



#' Plot frame coverage of 5'-site profiles
#'
#' This function plots the frame coverage of 5'-site profiles,
#' i.e. the fraction of reads (their 5' end) in frame 0, 1, and 2 (defined by start codon),
#' as the mean over all transcripts with its standard error
#' for individual read lengths as well as for all read lengths summarized.
#'
#' This plot is used in the RiboseQC report in section 4.2.1.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param analysis_frame_cutoff Object containing statistics on frame coverage
#' (per originating compartment and read length) generated by \code{RiboseQC_analysis}. \cr \cr
#' Example: \cr
#' \code{res_all$selection_cutoffs$analysis_frame_cutoff}
#'
#' @param comp String for originating compartment \cr \cr
#' Check for available originating compartments (with data) in the data set using:
#' \code{names(Filter(length, res_all$selection_cutoffs$analysis_frame_cutoff))}
#'
#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plot for the html report and returns \code{NULL}, invisibly.
#' With \code{output_rds_path}, it also saves the plot as the RDS object file \code{<output_rds_path><sample>_<comp>_4_selection5toP}.
#'
#' @seealso \code{\link{create_html_report}}
#' @export


plot_frame_dist_boxplot <- function(analysis_frame_cutoff, comp, sample="", output_rds_path="") {
    rl <- names(analysis_frame_cutoff[[comp]])
    rl_list <- NULL

    for (k in rl){
        data <- analysis_frame_cutoff[[comp]][[k]]$frames
        colnames(data) <- c("frame 0", "frame 1", "frame 2")
        rl_list[[k]] <- data
    }
    data <- suppressMessages(melt(rl_list))

    ## order read lengths
    data$L1 <- factor(data$L1, levels=sort(names(analysis_frame_cutoff[[comp]])))

    gg <- ggplot(data, aes(x = L1, y = value, color=Var2, group=Var2)) +
        stat_summary(fun.data = mean_se,size=1.5) +
        scale_x_discrete(name="read length") +
        scale_color_manual(values = alpha(c("red","forestgreen","blue"),alpha = .8)) +
        scale_y_continuous(name="read distribution between frames",limits=c(0,1)) +
        theme(legend.position = "bottom", panel.background=element_blank()) +
        labs(color = "frame") +
        guides(fill=guide_legend(title=""))


    # print to HTML file
    print(gg)

    cat("\n\n")  # necessary for RMarkdown rendering

    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_", comp, "_4_selection5toP")
        saveRDS(list(gg, c(7,5)), file = file_path,
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}


#' Get selected read lengths and cutoffs
#'
#' This function retrieves data on the read lengths with a cutoff, selected or not (per originating compartment),
#' e.g. their cutoff values, frame preference and codon gain. \cr \cr
#' This data is used in the RiboseQC report in section 4.2.2 (displayed as table).
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @return This function returns a list with a list for each sample, with a data.frame for each originating compartment:
#' the read lengths, their cutoffs, frame preference and codon gain.
#'
#' @seealso \code{\link{create_html_report}}
#' @export

get_rl_and_cutoffs <- function(rdata_list) {
    datasets <- NULL
    for (i in c(seq_along(rdata_list))) {
        res_all <- rdata_list[[names(rdata_list)[i]]]

        sample <- names(rdata_list)[[i]]
        comps <- names(res_all$selection_cutoffs$results_choice)

        for (comp in comps){
            d <- as.data.frame(res_all$selection_cutoffs$results_choice[[comp]]$data)[c(3, 1, 2, 4, 5)]
            names(d) <- c("read length", "cutoff", "frame preference (%)", "gain_codon?", "gain_new_codons?")
            datasets[[sample]][[comp]] <- d
        }
    }
    return(datasets)
}



#' Get default choice of read lengths
#'
#' This function returns the read lengths per originating compartment that the default method
#' (\code{readlength_choice_method = "max_coverage"} of \code{RiboseQC_analysis}) selects.
#' With this method, these read lengths build the basis for all P-site based calculations such as
#' metagene and codon usage analyses
#' (e.g. \code{plot_metagene_hm} and \code{plot_codon_usage_positional}); with another method,
#' \code{offsets_df} or \code{rescue_all_rls}, other read lengths do. \cr \cr
#' With the default method, the RiboseQC report shows these read lengths in section 4.2.3 (displayed as table);
#' with another method, it shows the read lengths of that method (\code{final_choice}).
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @return This function returns a list with a matrix for each sample: the compartments and their read lengths.
#'
#' @seealso \code{\link{create_html_report}}
#' @export

get_default_rl_selection <- function(rdata_list){
    datasets <- NULL
    for (i in seq_along(rdata_list)){
        res_all <- rdata_list[[names(rdata_list)[i]]]
        sample <- names(rdata_list)[[i]]
        comps <- names(res_all$selection_cutoffs$results_choice)

        rl_choices_default <- NULL
        for (j in comps){
            data <- res_all$selection_cutoffs$results_choice[[j]]$data
            data_default <- data[data[["max_coverage"]],][["read_length"]] # choose only read lengths with max_coverage=TRUE (default choice)
            rl_choices_default <- rbind(rl_choices_default, c(compartment=j,
                                                              `read lengths of choice (ordered by quality)`=paste(data_default, collapse = ", ")))
        }
        datasets[[sample]] <- rl_choices_default
    }
    return(datasets)
}



#' Get top 50 mapping positions
#'
#' This function retrieves data on the top 50 positions (nucleotide resolution) are listed
#' where most reads (their 5' end) map to, revealing possibly contaminating sequences.
#' The reads are counted for each read length, so a position can be listed once for each read length. \cr \cr
#' This data is used in the RiboseQC report in section 5 (displayed as table).
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @return This function returns data to be displayed as table in the html report: a list with a data.frame for each sample.
#'
#' @seealso \code{\link{create_html_report}}
#' @export

get_top50_mapping <- function(rdata_list) {
    datasets <- list()
    for (i in c(seq_along(rdata_list))){
        res_all <- rdata_list[[names(rdata_list)[i]]]
        data <- as.data.frame(res_all$sequence_analysis)
        data <- data[,c("score", "pct", "seqnames", "region", "gene_biotype", "gene", "seq")]
        names(data) <- c("read count", "library fraction (%)", "seqnames", "region", "gene biotype", "gene", "sequence")
        datasets[[names(rdata_list)[i]]] <- data
    }
    return(datasets)
}



#' Get top 50 abundant genes (CDS regions for protein coding genes)
#'
#' This function retrieves data on the top 50 abundant CDS regions for protein coding genes. \cr \cr
#' This data is used in the RiboseQC report in section 6 (displayed as table).
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @return This function returns data to be displayed as table in the html report: a list with a data.frame for each sample,
#' the 50 rows of \code{res_all$read_stats$counts_cds_genes} with the most reads (all rows, if it has fewer).
#'
#' @seealso \code{\link{create_html_report}}
#' @export

get_top50_cds_genes <- function(rdata_list) {
    datasets <- list()
    for (i in c(seq_along(rdata_list))){
        res_all <- rdata_list[[names(rdata_list)[i]]]
        rc_cds <- as.data.frame(res_all$read_stats$counts_cds_genes)
        rc_cds <- rc_cds[with(rc_cds, order(-reads, chromosome)), ][seq_len(min(50, nrow(rc_cds))),]
        datasets[[names(rdata_list)[i]]] <- rc_cds
    }
    return(datasets)
}

#' Get top 50 abundant genes (all genes)
#'
#' This function retrieves data on the top 100 abundant genes (the function name says 50). \cr \cr
#' The RiboseQC report shows the first 50 of these rows in section 6 (displayed as table).
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @return This function returns data to be displayed as table in the html report: a list with a data.frame for each sample,
#' the 100 rows of \code{res_all$read_stats$counts_all_genes} with the most reads (all rows, if it has fewer).
#'
#' @seealso \code{\link{create_html_report}}
#' @export

get_top50_all_genes <- function(rdata_list) {
    datasets <- list()
    for (i in c(seq_along(rdata_list))){
        res_all <- rdata_list[[names(rdata_list)[i]]]
        rc_all <- as.data.frame(res_all$read_stats$counts_all_genes)
        rc_all <- rc_all[with(rc_all, order(-reads, chromosome)), ][seq_len(min(100, nrow(rc_all))),]
        datasets[[names(rdata_list)[i]]] <- rc_all
    }
    return(datasets)
}


#' Get codon usage data (positional and bulk)
#'
#' This function processes codon usage data generated by \code{RiboseQC_analysis}, i.e.
#' codon usage
#' \itemize{
#' \item per codon position (i.e. positional codon usage) as well as
#' \item summed up over all positions* (i.e. bulk codon usage)
#' }
#' for a specific data type, originating compartment, and read length. The amino acids of the codons come from
#' \code{RiboseQC_analysis}: the genetic code of each circular chromosome, and code 1 for the other chromosomes
#' (see \code{prepare_annotation_files}). \cr \cr
#' This data is used as input in \code{plot_codon_usage_bulk}
#' to generate a bar plot, and in \code{plot_codon_usage_bulk_rmd} to iterativly generate plots for the RiboseQC report in section 8. \cr\cr
#' Please check \code{plot_codon_usage_positional} for positional codon usage (instead of bulk codon usage). \cr\cr
#' * Information on positions included in analysis: \cr
#' Based on P-site corrected reads, the codon usage within CDS regions for protein coding genes
#' is examplatory calculated
#' \itemize{
#' \item for the first 11 codons of the CDS (referred to as \emph{start}),
#' \item for 11 codons from the middle of the CDS (referred to as \emph{middle}), and
#' \item for the last 11 codons of the CDS (referred to as \emph{stop}).
#' }
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param data Object (list of lists) generated by \code{RiboseQC_analysis}: \code{res_all}.
#'
#' @param data_type String; select one of the following: \cr
#' \itemize{
#' \item \code{Codon_counts}: codon count in defined positions*, \cr
#' \item \code{P_sites_percodon}: P-sites count in defined positions*, or \cr
#' \item \code{P_sites_percodon_ratio}: ratio of P-sites counts to codon counts in defined positions*.
#' \item \code{E_sites_percodon}: E-sites count in defined positions*, or \cr
#' \item \code{E_sites_percodon_ratio}: ratio of E-sites counts to codon counts in defined positions*.
#' \item \code{A_sites_percodon}: A-sites count in defined positions*, or \cr
#' \item \code{A_sites_percodon_ratio}: ratio of A-sites counts to codon counts in defined positions*.
#' }
#'
#' @param comp String for originating compartment. \cr
#' Check for available originating compartments in the data set using: \cr
#' \code{names(res_all$profiles_P_sites$Codon_counts)}
#'
#' @param rl String for read length. \cr
#' Check for available read lengths in the data set using: \cr
#' \code{names(res_all$profiles_P_sites[[data_type]][[comp]])}
#'
#' @return This function returns a list (e.g. called codon_usage_data) with information on codon usage: \cr
#' \itemize{
#' \item \code{codon_usage_data$data_positional} contains the data on positional codon usage, with no scaling (\code{none}),
#' log2 scaling (\code{log2}) and z scoring (\code{zscore}),
#' \item \code{codon_usage_data$data_bulk} contains the data on bulk codon usage,
#' \item \code{codon_usage_data$data_type} saves the data type used (see parameter \code{data_type})
#' \item \code{codon_usage_data$comp} saves the originating compartment used (see parameter \code{comp}), and
#' \item \code{codon_usage_data$rl} saves the read length used (see parameter \code{rl}).
#' }
#'
#' @seealso \code{\link{create_html_report}}
#' @export

get_codon_usage_data <- function(data, data_type, comp, rl) {

    # test input parameters for existence
    if((data_type != "Codon_counts") & (data_type != "P_sites_percodon") & (data_type != "P_sites_percodon_ratio") & (data_type != "A_sites_percodon") & (data_type != "A_sites_percodon_ratio") & (data_type != "E_sites_percodon") & (data_type != "E_sites_percodon_ratio")){
        stop("Variable data_type takes only: \"Codon_counts\" or \"P_sites_percodon\" or \"P_sites_percodon_ratio\" or \"A_sites_percodon\" or \"A_sites_percodon_ratio\" \"E_sites_percodon\" or \"E_sites_percodon_ratio\"")
    }
    if(!(comp %in% names(data$profiles_P_sites[[data_type]]))) {
        stop("Variable comp takes only: ", paste(names(data$profiles_P_sites[[data_type]]), collapse=", "))
    }
    if(!(rl %in% names(data$profiles_P_sites[[data_type]][[comp]]))) {
        stop("Variable rl takes only: ", paste(names(data[["profiles_P_sites"]][[data_type]][[comp]]), collapse=", "))
    }

    # POSITIONAL CODON USAGE

    # get codon usage (per nucleotide position)
    d <- data[["profiles_P_sites"]][[data_type]][[comp]][[rl]]
    d <- as.data.frame(d)
    d[is.na(d)] <- 0 # set all NAs to zero
    names(d) <- seq(1,length(names(d)))

    # get row names (example row name: "GGG;G" corrsponding to codon;amino acid)
    row_names <- row.names(d)

    # get codons
    codons <- unlist(lapply(row_names, function(x) strsplit(x, ";")[[1]][1]))
    codons <- factor(codons, levels=rev(sort(codons)))

    # add amino acids
    gencode <- unlist(lapply(row_names, function(x) strsplit(x, ";")[[1]][2]))
    gencode <- factor(gencode, levels=c("M",
                                        "G", "A", "V", "I", "L", "F", "Y", "W",
                                        "C", "P",
                                        "S", "T", "N", "Q",
                                        "R", "H", "K",
                                        "D", "E",
                                        "*"))

    # add scaling (none, log2, zscore)
    # none
    d_none <- cbind(aa=gencode, codon=codons, d)
    # log2
    d_log2 <- cbind(aa=gencode, codon=codons, log2(d+1))
    # zscore
    d_zscore <- cbind(aa=gencode, codon=codons, data.frame(scale(d, center=TRUE, scale=TRUE)))
    colnames(d_zscore) <- c("aa", "codon", seq(1,length(names(d)))) # need to set position names again
    # save all scaling types as list
    data_positional <- list(none=d_none, log2=d_log2, zscore=d_zscore)

    # BULK CODON USAGE

    # summarize codon usage over all nucleotide positions
    d_bulk <- rowSums(d)
    d_bulk <- as.data.frame(d_bulk)
    d_bulk <- cbind(aa=gencode, codon=codons, d_bulk)

    # save codon usage bulk data and metadata information
    codon_usage_data <- list(data_positional=data_positional,
                             data_bulk=d_bulk,
                             data_type=data_type,
                             comp=comp,
                             rl=rl)

    return(codon_usage_data)
}


#' Plot positional codon usage heatmaps within the RMarkdown document
#'
#' This function generates iteratively all positional codon usage heatmaps;
#' iteration over originating compartments, read length, data type (codon count, and A-, P- and E-site counts and their ratios to the codon counts),
#' and scaling method (none, log2, zscale). \cr \cr
#' These plots are displayed in the RiboseQC report in section 7.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param data Object (list of lists) generated by \code{RiboseQC_analysis}: \code{res_all}.
#'
#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints iteratively all positional codon usage plots
#' for the html report and returns \code{NULL}. With \code{output_rds_path}, it also saves the same plots as RDS object files.
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' \dontrun{
#' # This draws all plots of section 7 of the report: 210 plots, about 7 min, on the bundled
#' # shoots sample (the root sample has no P-sites).
#' load(system.file("ext_data", "shoots_results_RiboseQC", package = "RiboseQC"))
#' plot_codon_usage_positional_rmd(res_all, "shoots")
#' }
#' @export

plot_codon_usage_positional_rmd <- function(data, sample="", output_rds_path="") {

    # different data types for codon usage available: codon count, read count, codon-read ratio
    data_types <- c("Codon counts"="Codon_counts",
                    "A-sites counts"="A_sites_percodon",
                    "P-sites counts"="P_sites_percodon",
                    "E-sites counts"="E_sites_percodon",
                    "A-sites per codon"="A_sites_percodon_ratio",
                    "P-sites per codon"="P_sites_percodon_ratio",
                    "E-sites per codon"="E_sites_percodon_ratio")

    # assumption that list structure is exactly the same for Codon_counts, P_sites_percodon, P_sites_percodon_ratio
    comps <- names(data[["profiles_P_sites"]][["Codon_counts"]])

    # per compartment
    for (comp in comps) {

        # cat commands necessary for RMarkdown rendering
        cat("\n\n")
        cat("### ", comp, " {.tabset} \n\n")
        cat("\n\n")
        rls <- as.vector(lapply(data[["profiles_P_sites"]][["Codon_counts"]], function(x) names(x))[[comp]])

        # per read length
        for (rl in rls) {
            cat("\n\n")
            cat("#### ", rl, " {.tabset} \n\n")
            cat("\n\n")

            # per codon count, read count, or codon-read ratio
            for(i in seq_along(data_types)){
                data_type <- data_types[i]
                cat("\n\n")
                cat("##### ", names(data_type), " {.tabset} \n\n")
                cat("\n\n")

                # get data
                codon_usage_data <- get_codon_usage_data(data, data_type, comp, rl)

                # per scaling method
                scals <- names(codon_usage_data$data_positional)
                for (scal in scals){
                    cat("\n\n")
                    cat("###### ", scal, " {.tabset} \n\n")
                    cat("\n\n")

                    # plot positional codon usage
                    plot_codon_usage_positional(codon_usage_data, scal, sample, output_rds_path)
                }
            }
        }
    }
}


#' Plot positional codon usage heatmap
#'
#' This function plots the codon usage per codon position* (positional codon usage) as heatmap
#' for a specific data type, originating compartment, read length, and scaling method. \cr \cr
#' * Information on positions included in analysis: \cr
#' Based on P-site corrected reads, the codon usage within CDS regions for protein coding genes
#' is examplatory calculated
#' \itemize{
#' \item for the first 11 codons of the CDS (referred to as \emph{start}),
#' \item for 11 codons from the middle of the CDS (referred to as \emph{middle}), and
#' \item for the last 11 codons of the CDS (referred to as \emph{stop}).
#' }
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param codon_usage_data List containing codon usage data and meta data, generated by \code{get_codon_usage_data}.
#'
#' @param scal Scaling method:
#' no scaling (\code{none}), log2 scaling (\code{log2}), or z scoring (\code{zscore}).
#'
#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plot for the html report and returns it, invisibly.
#' With \code{output_rds_path}, it also saves the plot as the RDS object file \code{<output_rds_path><sample>_<comp>_7_codonusage_positional_<data type>_<rl>_<scal>}.
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' # the P-sites of ChrC of the shoots sample: the root sample has no P-sites
#' load(system.file("ext_data", "shoots_results_RiboseQC", package = "RiboseQC"))
#' codon_usage_data <- get_codon_usage_data(res_all, "P_sites_percodon", "ChrC", "all")
#' plot_codon_usage_positional(codon_usage_data, "zscore", "shoots")
#' @export

plot_codon_usage_positional <- function(codon_usage_data, scal, sample="", output_rds_path="") {

    # get data and meta information
    d_positional <- codon_usage_data$data_positional[[scal]]
    data_type <- codon_usage_data$data_type
    comp <- codon_usage_data$comp
    rl <- codon_usage_data$rl

    # translate data_type into a string suitable for the y axis label
    data_types <- c("Codon counts"="Codon_counts",
                    "A-sites counts"="A_sites_percodon",
                    "P-sites counts"="P_sites_percodon",
                    "E-sites counts"="E_sites_percodon",
                    "A-sites per codon"="A_sites_percodon_ratio",
                    "P-sites per codon"="P_sites_percodon_ratio",
                    "E-sites per codon"="E_sites_percodon_ratio")
    data_type <- names(data_types[data_types==data_type])

    # since ATG contains much more data than the other codons,
    # we will plot two heatmaps: ATG and non-ATG

    # heatmap for start and stop codons
    d <- d_positional[d_positional$aa=="*" | d_positional$codon=="ATG", ]  # get start and stop codons
    d.m <- suppressMessages(melt(d))
    d.m$position <- factor(c(rep("start", nrow(d.m)/3),
                             rep("middle", nrow(d.m)/3),
                             rep("stop", nrow(d.m)/3)),
                           levels=c("start", "middle", "stop"))
    g1 <- ggplot(d.m, aes(x=variable, y=codon)) +
        geom_tile(aes(fill=value), colour="darkgrey") +
        facet_grid(aa~position, scales = "free", space = "free") +
        scale_fill_gradient(low="white", high="steelblue", na.value="gray80") +
        scale_y_discrete(expand=c(0, 0)) +
        theme(axis.text.y=element_text(size=18, family="mono"),
              axis.title.y=element_text(size=18),
              axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              axis.ticks=element_blank(),
              panel.background=element_blank(),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.title=element_text(size=18),
              legend.text=element_text(size=14),
              strip.text = element_text(size=18)) +
        labs(colour="r", y="", x="\nposition (codon resolution)") +
        geom_point(data=d.m, aes(size="NA"), shape=NA, colour="gray80") +
        guides(fill=guide_legend(title=paste0(data_type,"(scale: ", scal, ")"), nrow=1),
               size=guide_legend(title=""))
    # heatmap for NON start and stop codons
    d <- d_positional[!(d_positional$aa=="*" | d_positional$codon=="ATG"), ]  # remove start and stop codons
    d.m <- suppressMessages(melt(d))
    d.m$position <- factor(c(rep("start", nrow(d.m)/3),
                             rep("middle", nrow(d.m)/3),
                             rep("stop", nrow(d.m)/3)),
                           levels=c("start", "middle", "stop"))
    g2 <- ggplot(d.m, aes(x=variable, y=codon)) +
        geom_tile(aes(fill=value), colour="darkgrey") +
        facet_grid(aa~position, scales="free", space="free") +
        scale_fill_gradient(low="white", high="steelblue", na.value="gray80") +
        scale_y_discrete(expand=c(0, 0)) +
        theme(axis.text.x=element_text(size=18),
              axis.text.y=element_text(size=18, family="mono"),
              axis.title=element_text(size=18),
              axis.ticks = element_blank(),
              panel.background=element_blank(),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.title=element_text(size=18),
              legend.text=element_text(size=14,angle = 45),
              strip.text.y = element_text(size=18),
              strip.text.x = element_blank()) +
        labs(colour="r", y="", x="\nposition (codon resolution)") +
        geom_point(data=d.m, aes(size="NA"), shape=NA, colour="gray80") +
        labs(fill=paste0(data_type,"(scale: ", scal, ")")) + guides(size=guide_legend(title=""))

    # zscore: the plot with other colours below
    if(scal!="zscore"){gg <- arrangeGrob(g1, g2, ncol=1, heights = c(1, 8))}

    suppressMessages(if(scal=="zscore"){gg<-arrangeGrob(g1 + scale_fill_gradient2(low = "blue",mid = "white",high = "firebrick",na.value = "gray80",guide = "colourbar"),
                                       g2 + scale_fill_gradient2(low = "blue",mid = "white",high = "firebrick",na.value = "gray80",guide = "colourbar"),
                                       ncol=1, heights = c(1, 8))})

    cat("\n\n")  # necessary for RMarkdown rendering

    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_",
                           comp, "_",
                           "7_codonusage_positional_",
                           gsub(" ", "_", data_type), "_",
                           rl, "_",
                           scal)
        saveRDS(list(gg, c(12,16)), file = file_path,
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }

    # print to HTML file
    print(as_ggplot(gg))
}


#' Plot bulk codon usage bar plots within the RMarkdown document
#'
#' This function generates iteratively all bulk codon usage bar plots;
#' iteration over originating compartments, read length, and data type (codon count, and A-, P- and E-site counts and their ratios to the codon counts). \cr \cr
#' These plots are displayed in the RiboseQC report in section 8.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param data Object (list of lists) generated by \code{RiboseQC_analysis}: \code{res_all}.
#'
#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints iteratively all bulk codon usage plots
#' for the html report and returns \code{NULL}. With \code{output_rds_path}, it also saves the same plots as RDS object files.
#'
#' @seealso \code{\link{create_html_report}}, \code{\link{get_codon_usage_data}}, \code{\link{plot_codon_usage_bulk}}
#' @examples
#' # This draws all plots of section 8 of the report: 70, about 1 min, on the bundled shoots
#' # sample (the root sample has no P-sites).
#' load(system.file("ext_data", "shoots_results_RiboseQC", package = "RiboseQC"))
#' plot_codon_usage_bulk_rmd(res_all, "shoots")
#' @export

plot_codon_usage_bulk_rmd <- function(data, sample="", output_rds_path="") {

    # different data types for codon usage available: codon count, read count, codon-read ratio
    data_types <- c("Codon counts"="Codon_counts",
                    "A-sites counts"="A_sites_percodon",
                    "P-sites counts"="P_sites_percodon",
                    "E-sites counts"="E_sites_percodon",
                    "A-sites per codon"="A_sites_percodon_ratio",
                    "P-sites per codon"="P_sites_percodon_ratio",
                    "E-sites per codon"="E_sites_percodon_ratio")

    # assumption that list structure is exactly the same for Codon_counts, P_sites_percodon, P_sites_percodon_ratio
    comps <- names(data[["profiles_P_sites"]][["Codon_counts"]])

    for (comp in comps) {

        # cat commands necessary for RMarkdown rendering
        cat("\n\n")
        cat("### ", comp, " {.tabset} \n\n")
        cat("\n\n")
        rr <- as.vector(lapply(data[["profiles_P_sites"]][["Codon_counts"]], function(x) names(x))[[comp]])

        for (rl in rr) {
            cat("\n\n")
            cat("#### ", rl, " {.tabset} \n\n")
            cat("\n\n")

            for(i in seq_along(data_types)){
                data_type <- data_types[i]
                cat("\n\n")
                cat("##### ", names(data_type), " {.tabset} \n\n")

                # get data
                codon_usage_data <- get_codon_usage_data(data, data_type, comp, rl)

                # plot data
                plot_codon_usage_bulk(codon_usage_data, sample, output_rds_path)
            }
        }
    }
}


#' Plot bulk codon usage bar plots
#'
#' This function plots the codon usage summed up over all positions* (bulk codon usage) as bar plot
#' for a specific data type, originating compartment, and read length. \cr \cr
#' * Information on positions included in analysis: \cr
#' Based on P-site corrected reads, the codon usage within CDS regions for protein coding genes
#' is examplatory calculated
#' \itemize{
#' \item for the first 11 codons of the CDS (referred to as \emph{start}),
#' \item for 11 codons from the middle of the CDS (referred to as \emph{middle}), and
#' \item for the last 11 codons of the CDS (referred to as \emph{stop}).
#' }
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param codon_usage_data List containing codon usage data and meta data, generated by \code{get_codon_usage_data}.
#'
#' @param sample String; sample name (selected from the input names
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#'
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder, ending with "/" (the file name is appended to it).
#'
#' @return This function prints the plot for the html report and returns it, invisibly.
#' With \code{output_rds_path}, it also saves the plot as the RDS object file \code{<output_rds_path><sample>_<comp>_8_codonusage_bulk_<data type>_<rl>}.
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#' # the P-sites of ChrC of the shoots sample: the root sample has no P-sites
#' load(system.file("ext_data", "shoots_results_RiboseQC", package = "RiboseQC"))
#' codon_usage_data <- get_codon_usage_data(res_all, "P_sites_percodon", "ChrC", "all")
#' plot_codon_usage_bulk(codon_usage_data, "shoots")
#' @export

plot_codon_usage_bulk <- function(codon_usage_data, sample="", output_rds_path="") {

    d_bulk <- codon_usage_data$data_bulk
    d_bulk.m <- suppressMessages(melt(d_bulk))

    data_type <- codon_usage_data$data_type
    comp <- codon_usage_data$comp
    rl <- codon_usage_data$rl

    # translate data_type into a string suitable for the y axis label
    data_types <- c("Codon counts"="Codon_counts",
                    "A-sites counts"="A_sites_percodon",
                    "P-sites counts"="P_sites_percodon",
                    "E-sites counts"="E_sites_percodon",
                    "A-sites per codon"="A_sites_percodon_ratio",
                    "P-sites per codon"="P_sites_percodon_ratio",
                    "E-sites per codon"="E_sites_percodon_ratio")
    data_type <- names(data_types[data_types==data_type])

    # prepare subplots (grobs): g1 and g2
    g1 <- ggplot(d_bulk.m[match("ATG", d_bulk.m$codon), ], aes(x=codon, y=value)) +
        geom_bar(aes(fill=aa), stat="identity") +
        facet_grid(.~aa, scales="free", space="free") +
        labs(x="", y=data_type) +
        scale_y_continuous(expand=c(0, 0)) +  # remove y axis offset
        scale_fill_manual(values=c("gray")) + # same color for all amino acids
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=16, family="Courier"),
              axis.text.y=element_text(size=16),
              axis.title.x=element_blank(),
              axis.title.y=element_text(size=16),
              axis.ticks.x=element_blank(),
              panel.background=element_rect(),  # get panel background
              panel.border = element_blank(),  # remove borders from panel
              panel.grid.major = element_blank(),  # remove grid lines from panel
              panel.grid.minor = element_blank(),  # remove grid lines from panel
              strip.text=element_text(size=16),
              legend.position="none")

    g2 <- ggplot(d_bulk.m[-match("ATG", d_bulk.m$codon), ], aes(x=codon, y=value)) +
        geom_bar(aes(fill=aa), stat="identity") +
        facet_grid(.~aa, scales="free", space="free") +
        labs(x="", y=data_type) +
        scale_y_continuous(expand=c(0, 0)) +  # remove y axis offset
        scale_colour_hue() +  # hue colors
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=16, family="Courier"),
              axis.text.y=element_text(size=16),
              axis.title=element_blank(),
              axis.ticks.x=element_blank(),
              panel.background=element_rect(),  # get panel background
              panel.border = element_blank(),  # remove borders from panel
              panel.grid.major = element_blank(),  # remove grid lines from panel
              panel.grid.minor = element_blank(),  # remove grid lines from panel
              strip.text=element_text(size=16),
              legend.position="none")

    # arrange plots
    gg <- arrangeGrob(g1, g2, ncol=2, widths = c(1, 9))

    cat("\n\n")

    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_",
                           comp, "_",
                           "8_codonusage_bulk_",
                           gsub(" ", "_", data_type), "_",
                           rl)
        saveRDS(list(gg, c(12,8)), file = file_path,
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }

    # print to HTML file
    print(as_ggplot(gg))
}

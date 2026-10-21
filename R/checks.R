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

# The checks of the input, which run before the long steps (from ORFquant's R/checks.R).


#checks of the input, done before the long steps start, so that a wrong input stops the run with a clear message

check_files_exist<-function(files){
  files<-files[!is.na(files)]
  missing_files<-files[file.access(files,0)==-1]
  if(length(missing_files)>0){
    stop("The following files don't exist:\n",paste(missing_files,collapse="\n"),"\n",call. = FALSE)
  }
}

#' Checks of the input, for ORFquant
#'
#' \code{check_n_cores()} stops if \code{n_cores} is not one number, 1 or more.
#' \code{check_output_dir()} stops if the directory of \code{prefix} does not exist or cannot be
#' written to: the output files are written at the end of a run.
#' \code{\link{prepare_for_ORFquant}} uses them. RiboseQC exports them for ORFquant.
#' @param n_cores the number of cores
#' @param prefix the prefix of the output files, with their directory
#' @return \code{NULL}, invisibly, if the input is correct
#' @keywords internal
#' @name input_checks
NULL

#' @rdname input_checks
#' @export
check_n_cores<-function(n_cores){
  if(missing(n_cores) || !is.numeric(n_cores) || length(n_cores)!=1 || is.na(n_cores) || n_cores<1){
    stop("n_cores must be a number, 1 or more!",call. = FALSE)
  }
}

#the output files are written at the end of a run: their directory must exist and be writable
#' @rdname input_checks
#' @export
check_output_dir<-function(prefix){
  out_dir<-dirname(prefix)
  if(!dir.exists(out_dir) || file.access(out_dir,2)==-1){
    stop("Cannot write the output files, ",out_dir," is not a directory that exists and can be written to (",prefix,")",call. = FALSE)
  }
}

#the chromosomes of the BAM file must be chromosomes of the annotation, with the same lengths: else no read is on the
#annotation, or the first chunk of reads stops the run with "incompatible seqlengths"
check_bam_chromosomes<-function(bam_seqinfo,annotation_seqinfo,bam_file,annotation_file){
  if(!any(seqlevels(bam_seqinfo)%in%seqlevels(annotation_seqinfo))){
    stop("No chromosome of the BAM file is in the annotation: the BAM file ",bam_file," has e.g. ",paste(head(seqlevels(bam_seqinfo),3),collapse=", "),
         " and the annotation ",annotation_file," has e.g. ",paste(head(seqlevels(annotation_seqinfo),3),collapse=", "),
         ". The reads must be aligned to the genome of the annotation (not to its transcripts), with the same chromosome names",call. = FALSE)
  }
  common_chroms<-intersect(seqlevels(bam_seqinfo),seqlevels(annotation_seqinfo))
  diff_len<-common_chroms[which(seqlengths(bam_seqinfo)[common_chroms]!=seqlengths(annotation_seqinfo)[common_chroms])]
  if(length(diff_len)>0){
    stop("The BAM file ",bam_file," is not on the genome of the annotation ",annotation_file,": ",length(diff_len),
         " chromosomes have other lengths, e.g. ",diff_len[1]," has ",seqlengths(bam_seqinfo)[[diff_len[1]]]," nt in the BAM file and ",
         seqlengths(annotation_seqinfo)[[diff_len[1]]]," nt in the annotation",call. = FALSE)
  }
}

#the first exon line of the GTF file, which must have transcript_id and gene_id, tells the file apart from an empty
#one, or one without exon lines; this stops reading at that line
check_gtf<-function(gtf_file){
  con<-file(gtf_file,"r")
  on.exit(close(con))
  repeat{
    lines<-readLines(con,n=10000)
    if(length(lines)==0){break}
    exon_line<-grep("^[^\t]*\t[^\t]*\texon\t",lines,value = TRUE)
    if(length(exon_line)>0){
      if(!grepl("(^|[;\t ])transcript_id[ =]",exon_line[1]) || !grepl("(^|[;\t ])gene_id[ =]",exon_line[1])){
        stop("An exon line of the GTF file ",gtf_file," has no transcript_id or gene_id: the annotation needs both on every line",call. = FALSE)
      }
      return(invisible(TRUE))
    }
  }
  stop("The GTF file ",gtf_file," is empty or has no exon lines: the annotation needs the exons and the coding sequences (CDS) of the transcripts",call. = FALSE)
}

#a transcript gets the value of an attribute from its first line with one, and RiboseQC_analysis() gives a gene the
#name and biotype of its first transcript in trann: when the lines do not agree, their order in the GTF file chooses
#the value. gtf_ids: the attributes of the lines; trann: one row for each transcript
check_gtf_values<-function(gtf_ids,trann){
  #the ids (column by) with more than one value in a column of cols, and the column
  disagree<-function(x,by,cols,with_na){
    res<-lapply(cols,function(cl){
      v<-unique(x[!is.na(x[,by]) & (with_na | !is.na(x[,cl])),c(by,cl)])
      sort(unique(v[duplicated(v[,by]),by]))
    })
    data.frame(id=unlist(res),attribute=rep(cols,lengths(res)),stringsAsFactors=FALSE)
  }
  examples<-function(d){paste(head(paste0(d$id," (",d$attribute,")"),3),collapse=", ")}
  txs<-disagree(gtf_ids,"transcript_id",setdiff(names(gtf_ids),"transcript_id"),with_na=FALSE)
  genes<-disagree(trann,"gene_id",c("gene_name","gene_biotype"),with_na=TRUE)
  msg<-NULL
  if(nrow(txs)>0){
    n<-length(unique(txs$id))
    msg<-c(msg,paste0("the lines of ",n," ",ngettext(n,"transcript","transcripts")," give different values (e.g. ",examples(txs),
                      "), and a transcript gets the value of its first line with one"))
  }
  if(nrow(genes)>0){
    n<-length(unique(genes$id))
    msg<-c(msg,paste0("the transcripts of ",n," ",ngettext(n,"gene","genes")," give different gene names or biotypes (e.g. ",examples(genes),
                      "), and RiboseQC_analysis() gives a gene the values of its first transcript in the file"))
  }
  if(length(msg)>0){
    warning("The lines of the GTF file do not agree, so their order chooses some values: ",paste(msg,collapse="; "),
            " (see ?prepare_annotation_files)",call. = FALSE)
  }
}

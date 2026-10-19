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

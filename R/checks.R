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

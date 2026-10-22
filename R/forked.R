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

# Forked processes for steps that only read their data (from ORFquant's R/helpers.R).


#' Forked processes, for ORFquant
#'
#' \code{mc_job(expr, n_cores)} evaluates \code{expr} in a forked process
#' (\code{\link[parallel]{mcparallel}}) if \code{n_cores} is more than 1, and in this process if not.
#' \code{mc_value(job)} waits for the job and gives the value of \code{expr}. With
#' \code{wait = FALSE}, it gives \code{NULL} if the job is not done yet.
#' \code{mc_lapply(fs, n_cores)} is \code{lapply(fs, function(f) f())}, with each \code{f()} in a
#' forked process, at most \code{n_cores} at a time (\code{\link[parallel]{mclapply}}).
#' The warnings of a forked process are given again in this process, and its error stops here.
#' \code{\link{prepare_for_ORFquant}} and \code{\link{RiboseQC_analysis}} use \code{mc_job()} and \code{mc_value()}.
#' RiboseQC exports them for ORFquant.
#' @param expr the expression to evaluate
#' @param n_cores the number of cores
#' @param job a value of \code{mc_job()}
#' @param wait \code{FALSE} to not wait for the job
#' @param fs a list of functions without arguments
#' @return \code{mc_job()}: a job for \code{mc_value()}. \code{mc_value()}: the value of
#' \code{expr}. \code{mc_lapply()}: the list of the values of the functions.
#' @keywords internal
#' @name forked_processes
NULL

#' @rdname forked_processes
#' @export
mc_job<-function(expr,n_cores){
  if(n_cores>1){
    return(parallel::mcparallel(mc_capture(function(){expr})))
  }
  list(value=expr)
}
#' @rdname forked_processes
#' @export
mc_value<-function(job,wait=TRUE){
  if(!inherits(job,"parallelJob")){
    return(job$value)
  }
  res<-parallel::mccollect(job,wait=wait)
  if(is.null(res)){
    return(NULL)
  }
  mc_release(res[[1]])
}
#' @rdname forked_processes
#' @export
mc_lapply<-function(fs,n_cores){
  if(n_cores>1){
    return(lapply(parallel::mclapply(fs,mc_capture,mc.cores=n_cores,mc.preschedule=FALSE),mc_release))
  }
  lapply(fs,function(f){f()})
}
mc_capture<-function(f){
  warns<-list()
  value<-withCallingHandlers(f(),warning=function(w){
    warns[[length(warns)+1]]<<-w
    invokeRestart("muffleWarning")
  })
  list(value=value,warnings=warns)
}
mc_release<-function(res){
  if(inherits(res,"try-error")){
    stop(attr(res,"condition"))
  }
  if(is.null(res)){
    stop("A forked process ended without a result (out of memory?)")
  }
  for(w in res$warnings){
    warning(w)
  }
  res$value
}

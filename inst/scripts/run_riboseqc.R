## Run the two steps of RiboseQC from the command line, each in a new R
## process: prepare_annotation_files() and RiboseQC_analysis(). Run it with
## Rscript; --help lists the options. The installed copy is at
## system.file("scripts", "run_riboseqc.R", package = "RiboseQC").

usage <- "Usage: Rscript run_riboseqc.R (--gtf FILE --fasta FILE | --annotation FILE)
                             --bam FILE[,FILE...] --outdir DIR [options]

Runs prepare_annotation_files() and RiboseQC_analysis(), each in a new R
process, and writes all output files to DIR. The other parameters of the two
functions have their default values.

  --gtf FILE             annotation in GTF format (may be gzipped)
  --fasta FILE           genome sequence in FASTA format (may be bgzipped); the
                         annotation refers to it by its path, so keep it there
  --annotation FILE      the _Rannot file of an earlier run (in its
                         DIR/annotation), instead of --gtf and --fasta
  --bam LIST             Ribo-seq reads aligned to the genome with a spliced
                         aligner (e.g. STAR), sorted and indexed: a BAM file,
                         or a comma-separated list of BAM files
  --outdir DIR           output directory, created if it doesn't exist
  --sample-names LIST    comma-separated names of the samples, one for each BAM
                         file: the prefixes of their output files and their
                         names in the report (default: the BAM files' names
                         without .bam)
  --stranded VALUE       TRUE (default), FALSE, or inverse for libraries with
                         the reads on the strand opposite to the RNA; the same
                         for all BAM files
  --readlength-method M  how to select the read lengths for the P-sites:
                         max_coverage (default), max_inframe or all
  --cores N              number of cores for RiboseQC_analysis(), each for one
                         chunk of alignments in a forked process (default: 1)
  --no-report            don't create the HTML report
  --extended-report      create the extended HTML report, with the codon
                         occupancy of each read length
  -h, --help             print this help and exit

With --gtf and --fasta, the annotation is written to DIR/annotation; give its
_Rannot file to --annotation to reuse it for other samples. The output files of
a sample start with DIR/<its name>_, and the report of all samples is
DIR/<first name>_RiboseQC_report.html. At the end, the output files are listed,
one per line: their kind, a tab and their path.
"

args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("-h", "--help"))) {
  cat(usage)
  quit(status = 0)
}
if (length(args) == 0) {
  cat(usage)
  quit(status = 1)
}

flags <- c("no-report", "extended-report")
opts <- c(gtf = NA, fasta = NA, annotation = NA, bam = NA, outdir = NA, "sample-names" = NA,
          stranded = "TRUE", "readlength-method" = "max_coverage", cores = "1",
          "no-report" = "FALSE", "extended-report" = "FALSE")
i <- 1
while (i <= length(args)) {
  key <- sub("=.*", "", sub("^--", "", args[i]))
  if (!startsWith(args[i], "--") || !key %in% names(opts)) {
    stop("unknown option ", args[i], "; see --help")
  }
  if (key %in% flags) {
    if (grepl("=", args[i])) {
      stop("--", key, " takes no value")
    }
    opts[[key]] <- "TRUE"
  } else if (grepl("=", args[i])) {
    opts[[key]] <- sub("^[^=]*=", "", args[i])
  } else {
    if (i == length(args)) {
      stop("no value for ", args[i])
    }
    i <- i + 1
    opts[[key]] <- args[i]
  }
  i <- i + 1
}

from_gtf <- !is.na(opts[["gtf"]]) || !is.na(opts[["fasta"]])
if (from_gtf == !is.na(opts[["annotation"]])) {
  stop("give either --gtf and --fasta, or --annotation")
}
needed <- c(if (from_gtf) c("gtf", "fasta") else "annotation", "bam", "outdir")
if (any(is.na(opts[needed]))) {
  stop("missing ", paste0("--", needed[is.na(opts[needed])], collapse = ", "))
}
for (f in setdiff(needed, c("bam", "outdir"))) {
  if (!file.exists(opts[[f]])) {
    stop("--", f, ": file not found: ", opts[[f]])
  }
  opts[[f]] <- normalizePath(opts[[f]])
}
split_list <- function(x) {
  trimws(strsplit(x, ",")[[1]])
}
bam_files <- split_list(opts[["bam"]])
for (f in bam_files) {
  if (!file.exists(f)) {
    stop("--bam: file not found: ", f)
  }
}
bam_files <- normalizePath(bam_files)
if (is.na(opts[["sample-names"]])) {
  samples <- sub("\\.bam$", "", basename(bam_files))
} else {
  samples <- split_list(opts[["sample-names"]])
  if (length(samples) != length(bam_files)) {
    stop("--sample-names must have one name for each of the ", length(bam_files),
         " BAM files, not ", length(samples))
  }
}
if (anyDuplicated(samples)) {
  stop("two BAM files have the sample name ", samples[duplicated(samples)][1],
       "; give other names with --sample-names")
}
if (!opts[["stranded"]] %in% c("TRUE", "FALSE", "inverse")) {
  stop("--stranded must be TRUE, FALSE or inverse, not ", opts[["stranded"]])
}
if (!opts[["readlength-method"]] %in% c("max_coverage", "max_inframe", "all")) {
  stop("--readlength-method must be max_coverage, max_inframe or all, not ",
       opts[["readlength-method"]])
}
if (!grepl("^[1-9][0-9]*$", opts[["cores"]])) {
  stop("--cores must be a whole number above 0, not ", opts[["cores"]])
}
create_report <- opts[["no-report"]] == "FALSE"
if (!create_report && opts[["extended-report"]] == "TRUE") {
  stop("give either --no-report or --extended-report")
}
if (!nzchar(system.file(package = "RiboseQC"))) {
  stop("RiboseQC is not installed in this R")
}
## Each step runs in a new R process, so that RiboseQC_analysis() and its
## forked processes don't keep the memory of the annotation step. library()
## first: loaded by RiboseQC::f(), RiboseQC makes R warn "stack imbalance".
run_step <- function(call) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("-e", shQuote("suppressPackageStartupMessages(library(RiboseQC))"),
                      "-e", shQuote(deparse1(call))))
  if (status != 0) {
    stop(deparse(call[[1]]), "() failed", call. = FALSE)
  }
}

dir.create(opts[["outdir"]], recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(opts[["outdir"]])
prefixes <- file.path(outdir, samples)
report_file <- paste0(prefixes[1], "_RiboseQC_report.html")

if (from_gtf) {
  annotation_dir <- file.path(outdir, "annotation")
  run_step(bquote(RiboseQC::prepare_annotation_files(annotation_directory = .(annotation_dir),
                                                     gtf_file = .(opts[["gtf"]]),
                                                     genome_seq = .(opts[["fasta"]]))))
  annotation_file <- file.path(annotation_dir, paste0(basename(opts[["gtf"]]), "_Rannot"))
} else {
  annotation_file <- opts[["annotation"]]
}

stranded <- if (opts[["stranded"]] == "inverse") "inverse" else as.logical(opts[["stranded"]])
run_step(bquote(RiboseQC::RiboseQC_analysis(annotation_file = .(annotation_file),
                                            bam_files = .(bam_files),
                                            readlength_choice_method = .(opts[["readlength-method"]]),
                                            stranded = .(stranded),
                                            dest_names = .(prefixes),
                                            sample_names = .(samples),
                                            create_report = .(create_report),
                                            report_file = .(report_file),
                                            extended_report = .(opts[["extended-report"]] == "TRUE"),
                                            n_cores = .(as.integer(opts[["cores"]])))))

sample_files <- c(results_all = "_results_RiboseQC_all", results = "_results_RiboseQC",
                  offsets = "_P_sites_offsets.tsv", P_sites_calcs = "_P_sites_calcs",
                  for_ORFquant = "_for_ORFquant", junctions = "_junctions",
                  P_sites_plus = "_P_sites_plus.bedgraph", P_sites_minus = "_P_sites_minus.bedgraph",
                  P_sites_uniq_plus = "_P_sites_uniq_plus.bedgraph",
                  P_sites_uniq_minus = "_P_sites_uniq_minus.bedgraph",
                  coverage_plus = "_coverage_plus.bedgraph", coverage_minus = "_coverage_minus.bedgraph",
                  coverage_uniq_plus = "_coverage_uniq_plus.bedgraph",
                  coverage_uniq_minus = "_coverage_uniq_minus.bedgraph",
                  summary = "_summary.tsv", read_lengths = "_read_lengths.tsv", regions = "_regions.tsv",
                  genes = "_genes.tsv", P_sites_calcs_tsv = "_P_sites_calcs.tsv",
                  top50_positions = "_top50_positions.tsv")
outputs <- c(annotation = annotation_file,
             unlist(lapply(prefixes, function(p) setNames(paste0(p, sample_files), names(sample_files)))),
             report = report_file,
             report_plots = paste0(report_file, "_plots"))
outputs <- outputs[file.exists(outputs)]
cat("\nRiboseQC output files:\n")
cat(paste(names(outputs), outputs, sep = "\t"), sep = "\n")

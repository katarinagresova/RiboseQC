## How inst/ext_data/root_results_RiboseQC was made, for the examples of the
## help pages: RiboseQC_analysis() with its default parameters on the bundled
## root BAM file, without the HTML report. RiboseQC_analysis() writes the file
## as <dest_name>_results_RiboseQC; the file holds one object, res_all.
##
## Run from the package's root directory, with RiboseQC installed.

library(RiboseQC)
ext_data <- system.file("ext_data", package = "RiboseQC")
out_dir <- tempfile("root_results_")
dir.create(out_dir)

annotation_file <- prepare_annotation_files(annotation_directory = out_dir,
                                            gtf_file = file.path(ext_data, "simp_arab.gtf.gz"),
                                            genome_seq = file.path(ext_data, "simp_arab.fa.gz"))
RiboseQC_analysis(annotation_file = annotation_file,
                  bam_files = file.path(ext_data, "simp_arab_root.bam"),
                  dest_names = file.path(out_dir, "root"),
                  create_report = FALSE)
file.copy(file.path(out_dir, "root_results_RiboseQC"), "inst/ext_data", overwrite = TRUE)

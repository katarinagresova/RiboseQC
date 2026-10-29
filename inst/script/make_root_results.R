## How inst/ext_data/root_results_RiboseQC, inst/ext_data/shoots_results_RiboseQC and
## inst/ext_data/shoots_five_prime_subcodon were made, for the examples of the help pages:
## RiboseQC_analysis() with its default
## parameters on the bundled root and shoots BAM files, without the HTML report.
## RiboseQC_analysis() writes each file as <dest_name>_results_RiboseQC; each file holds
## one object, res_all. The root sample has no P-sites (no read length has an offset);
## the shoots sample has the P-sites of ChrC.
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
                  bam_files = file.path(ext_data, c("simp_arab_root.bam", "simp_arab_shoots.bam")),
                  dest_names = file.path(out_dir, c("root", "shoots")),
                  create_report = FALSE)
file.copy(file.path(out_dir, c("root_results_RiboseQC", "shoots_results_RiboseQC")), "inst/ext_data",
          overwrite = TRUE)
## inst/ext_data/shoots_five_prime_subcodon, for the examples of calc_cutoffs_from_profiles()
## and choose_readlengths(): the 5' profiles of each transcript of the shoots sample. They
## are only in the _results_RiboseQC_all file (res_all$profiles_fivepr$five_prime_subcodon).
## The file holds one object, five_prime_subcodon.
load(file.path(out_dir, "shoots_results_RiboseQC_all"))
five_prime_subcodon <- res_all$profiles_fivepr$five_prime_subcodon
save(five_prime_subcodon, file = "inst/ext_data/shoots_five_prime_subcodon")

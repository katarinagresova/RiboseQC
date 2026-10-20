# RiboseQC
A comprehensive analysis tool for Ribo-seq and small RNA-seq data


**Ribo-seQC** (*RiboseQC*) is an R package that performs quality control analysis of small RNA-seq data, with a focus on Ribo-seq and related techniques. Thanks to syntax and functions present in *Bioconductor* packages like *GenomicFeatures*, *rtracklayer* or *BSgenome*, this package can perform comprehensive analyses on a variety of genomic regions. In addition, Ribo-seQC automatically generates an HTML report for the analyzed samples, allowing for quick and interactive comparison of multiple samples at once.

This tool focuses on the analysis of different read lengths, taking into account the genomic regions they map to (e.g. coding sequence, UTRs, non-coding RNAs, mitochondria or chloroplasts). It calculates the P-site offset of each read length from the 5' ends of the reads around start codons, selects the read lengths with in-frame signal, and gives their P-sites. Other useful features, such as codon usage and the analysis of the top mapping positions, are available in the Ribo-seQC package.

This repository is a fork of [ohlerlab/RiboseQC](https://github.com/ohlerlab/RiboseQC) that installs and runs on current Bioconductor. [NEWS.md](NEWS.md) lists what changed, including fixes that change some results.

More details can be found in our manuscript:

### Ribo-seQC: comprehensive analysis of cytoplasmic and organellar ribosome profiling data ###

*Lorenzo Calviello, Dominique Sydow, Dermot Harnett, Uwe Ohler*

**bioRxiv** (2019), doi: https://doi.org/10.1101/601468


## Installation

*RiboseQC* needs R 4.5 or later and Bioconductor 3.22 or later.

With *BiocManager*:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("katarinagresova/RiboseQC")
```

The HTML report also needs [pandoc](https://pandoc.org). On Linux, the packages are usually built from source, which needs the development files of system libraries such as libcurl (`libcurl-devel` on RHEL or Fedora, `libcurl4-openssl-dev` on Debian or Ubuntu).

Or with conda, using [environment.yml](environment.yml) from this repository (R 4.5, Bioconductor 3.22, all dependencies and pandoc):

```bash
conda env create -f environment.yml
conda activate riboseqc
Rscript -e 'remotes::install_github("katarinagresova/RiboseQC", upgrade = "never")'
```

To also install the vignette, add `build_vignettes = TRUE` to either install command (this needs pandoc), and read it with `vignette("RiboseQC")`.


## Usage

Two steps are required to use *RiboseQC*: `prepare_annotation_files()` reads a GTF file and the genome sequence, once for each annotation, and `RiboseQC_analysis()` analyzes one or more BAM files with this annotation and creates the HTML report. This example runs both steps on the data included in the package: Ribo-seq reads from *Arabidopsis thaliana* roots and shoots (Hsu *et al*, PNAS 2016) on 309 genes of the nuclear, mitochondrial (`ChrM`) and chloroplast (`ChrC`) chromosomes.

```r
library(RiboseQC)
ext_data <- system.file("ext_data", package = "RiboseQC")
dir.create("example")

# The bundled genome sequence is compressed with gzip, which Rsamtools can't
# read: write an uncompressed copy
genome_fasta <- file.path(getwd(), "example", "simp_arab.fa")
Biostrings::writeXStringSet(Biostrings::readDNAStringSet(file.path(ext_data, "simp_arab.fa.gz")),
                            genome_fasta)

# 1. Once for each annotation: writes example/annotation/simp_arab.gtf.gz_Rannot
annotation_file <- prepare_annotation_files(annotation_directory = "example/annotation",
                                            gtf_file = file.path(ext_data, "simp_arab.gtf.gz"),
                                            genome_seq = genome_fasta)

# 2. Analyze the two BAM files: writes example/root_* and example/shoots_*,
#    and one HTML report for both samples, example/report.html
RiboseQC_analysis(annotation_file = annotation_file,
                  bam_files = file.path(ext_data, c("simp_arab_root.bam", "simp_arab_shoots.bam")),
                  dest_names = c("example/root", "example/shoots"),
                  sample_names = c("root", "shoots"),
                  report_file = "example/report.html")
```

The report can be opened using a browser like Chrome, Firefox or others. It includes all plots, so it is large: 60 MB for this example. The plots are also in `example/report.html_plots/`, as PDF files (`pdf/`) and as R objects (`rds/`: `readRDS(file)[[1]]` is a ggplot, or a gtable to draw with `grid::grid.draw()`). `create_html_report()` creates a new report from the `*_results_RiboseQC` files of one or more analyses.

The [vignette](vignettes/RiboseQC.Rmd) runs the same analysis without the report, and shows some of the results: the reads of each compartment, the read lengths, the P-site offsets and the frame of the P-sites. To read it with its results, without running it, open [vignettes/rendered/RiboseQC.md](vignettes/rendered/RiboseQC.md).

For your own data, give your GTF, FASTA and BAM files, and a `dest_names` prefix for each BAM file. Without `dest_names`, the output files go next to the BAM files. `?RiboseQC_analysis` describes all parameters, for example `stranded`, `readlength_choice_method` and `n_cores`.


### Input files

- **GTF**: `exon` and `CDS` lines, each with `transcript_id` and `gene_id`. `gene_name`, `gene_biotype` (or `gene_type`) and `transcript_biotype` (or `transcript_type`) are used when present.
- **Genome sequence**: an uncompressed FASTA file, with the same chromosome names as the GTF and the BAM files. Give its absolute path: the `_Rannot` file refers to it by the path that you give, and `RiboseQC_analysis()` reads it from there. Chromosomes named as in `circ_chroms` (by default `chrM`, `MT`, `ChrM`, `ChrC`, `Pltd`, `Pt` and other names of mitochondria and plastids) are treated as circular, and each one is analyzed as a separate compartment. `prepare_annotation_files()` no longer forges a *BSgenome* package from a 2bit file, which failed with current Bioconductor (ohlerlab/RiboseQC#21, #23); `twoBitToFa` of the UCSC tools (http://hgdownload.soe.ucsc.edu/admin/exe/) makes a FASTA file from a 2bit file.
- **BAM**: Ribo-seq reads aligned to the genome (not to the transcriptome) with a spliced aligner (e.g. STAR), sorted by coordinate and indexed with `samtools index`. Reads flagged as duplicates, secondary alignments and alignments with an insertion or a deletion (`I` or `D` in the CIGAR) are not used. Soft-clipped bases do not count in the read length, so remove adapters and UMIs before the alignment. Reads with a mapping quality above 50 count as uniquely mapping (STAR gives them 255).


### Output files

For each BAM file, with its entry of `dest_names` as the prefix:

- `<dest_name>_results_RiboseQC_all`: all results, in the list `res_all` (to read with `load()`): the read length distribution of each compartment, the reads in each genomic region and biotype, the counts, RPKM and TPM values of each gene, the profiles of 5' ends and P-sites, the offsets, the P-sites, the codon usage and the 50 top mapping positions. The *Details* section of `?RiboseQC_analysis` describes each element. With `write_tmp_files = FALSE`, this file is not written.
- `<dest_name>_results_RiboseQC`: a smaller copy of `res_all` for the report (`create_html_report()`): without the P-sites and the positions of the reads, and with the profiles summed over the transcripts.
- `<dest_name>_P_sites_offsets.tsv`: the offsets used for the P-sites, as a tab-separated table with the columns `read_length`, `cutoff` (the distance of the P-site from the read's 5' end) and `compartment` (`nucl` for all chromosomes that are not circular, or the name of a circular one).
- `<dest_name>_P_sites_calcs`: the read lengths of the analysis in each compartment (`comp`), with their offset (`cutoff`), frame preference and percentage of the compartment's reads (`pct_map`), and for each method of `readlength_choice_method` (`all`, `max_inframe`, `max_coverage`) if it selects them.
- `<dest_name>_P_sites_plus.bedgraph` and `<dest_name>_P_sites_minus.bedgraph`: the number of P-sites at each position, from all reads; `<dest_name>_P_sites_uniq_*.bedgraph`: the same from the uniquely mapping reads.
- `<dest_name>_coverage_plus.bedgraph`, `<dest_name>_coverage_minus.bedgraph` and `<dest_name>_coverage_uniq_*.bedgraph`: the coverage of the reads. With `normalize_cov = TRUE` (the default), the coverage is scaled so that the two strands sum to about one million.
- `<dest_name>_junctions`: the reads on the annotated splice junctions.
- `<dest_name>_for_ORFquant`: the P-sites and the junction reads for *ORFquant* (see below).

The HTML report is `report_file` (by default, the first BAM file followed by `_RiboseQC_report.html`), with its plots in `<report_file>_plots/`.


### ORFquant

[ORFquant](https://github.com/katarinagresova/ORFquant) detects and quantifies the translation of ORFs, and uses *RiboseQC* for its annotation step: the `_Rannot` file of `prepare_annotation_files()` is also the annotation of `run_ORFquant()`. `run_ORFquant()` reads the P-sites from `<dest_name>_for_ORFquant`. `prepare_for_ORFquant()` creates the same file from a BAM file and the offsets in `<dest_name>_P_sites_offsets.tsv`.


### Notes

- `library(RiboseQC)` attaches only *GenomicRanges* (and the packages it attaches). Scripts that use, for example, *Biostrings*, *rtracklayer* or *ggplot2* need their own `library()` calls.
- `n_cores` above 1 uses forked processes, which need Linux or macOS.
- `prepare_annotation_files()` needs temporary space for the *TxDb*. If the temporary directory is too small, for example on some compute clusters, it stops with `database or disk is full`; set `TMPDIR` to a directory with more space before starting R.


For problems with this fork, please open an issue at https://github.com/katarinagresova/RiboseQC/issues.

For any question, please email:

calviello.l.bio@gmail.com, dominique.sydow@posteo.de (analysis and data visualization), Dermot.Harnett@mdc-berlin.de (package maintainer), Uwe.ohler@mdc-berlin.de (project supervisor), katarina.gresova@mdc-berlin.de (this fork).


Enjoy!

# Example data

A small data set for the vignette, the examples of the help pages and the
tests. It covers 309 protein-coding genes of *Arabidopsis thaliana*: 99 on the
five nuclear chromosomes, 88 on the chloroplast genome (`ChrC`) and 122 on the
mitochondrial genome (`ChrM`). The GTF, FASTA and BAM files were made in 2019.

| File | Content | Source |
|---|---|---|
| `simp_arab.gtf.gz` | The `exon` and `CDS` lines of the 309 genes (410 transcripts). They come from a version of the Araport11 annotation, with the non-coding RNAs of TAIR10, that the authors of RiboseQC changed (the source column is `Araport_11_mod_Lorenzo`). The changes are not documented. | Araport11 [2] |
| `simp_arab.fa.gz`, `.fai`, `.gzi` | The TAIR10 genome sequence (`Chr1` to `Chr5`, `ChrC` and `ChrM`). Only the exons of the 309 genes with 100 bp on each side, and 4 regions on `Chr4` and `Chr5` (5 kb), keep their bases. All other bases (99.6% of the genome) are N. Compressed with `bgzip` and indexed, so that *Rsamtools* can read it. | TAIR10 [3, 4] |
| `simp_arab_root.bam`, `simp_arab_shoots.bam` and their `.bai` indexes | Ribo-seq alignments of a root sample (26,702 alignments) and a shoot sample (53,142). Most of them are on the 309 genes. The others (1,702 and 22,585) are on nuclear regions more than 50 kb from these genes, where the genome sequence is N. RiboseQC counts them as intergenic. Each alignment keeps its position, flag, mapping quality and CIGAR, but not its read name, sequence, base qualities or tags. The mapping qualities (255, 3, 1 and 0) are those of STAR. | Hsu et al. [1] |
| `root_results_RiboseQC`, `shoots_results_RiboseQC` | The `_results_RiboseQC` files of `RiboseQC_analysis()` with its default parameters on `simp_arab_root.bam` and `simp_arab_shoots.bam`. Each holds one object, `res_all`. The examples of the help pages use them. The root sample has no P-sites, because no read length has an offset; the shoots sample has the P-sites of `ChrC`. `inst/script/make_root_results.R` makes them. | RiboseQC |
| `shoots_five_prime_subcodon` | The 5' profiles of each transcript of the shoots sample, from the same run: `res_all$profiles_fivepr$five_prime_subcodon` of its `_results_RiboseQC_all` file. It holds one object, `five_prime_subcodon`. The examples of `calc_cutoffs_from_profiles()` and `choose_readlengths()` use it. `inst/script/make_root_results.R` makes it. | RiboseQC |

The reads are from GEO series GSE81332 (SRA SRP074840, BioProject
PRJNA321304), "Ribosome profiling in Arabidopsis root and shoot" [1]. It has
three Ribo-seq replicates of roots (GSM2150567, GSM2150568, GSM2150569) and
three of shoots (GSM2150570, GSM2150571, GSM2150572). The files do not record
which replicate each BAM file comes from, or if the replicates were pooled.

`inst/script/creating_test_data.R` shows how the GTF, FASTA and BAM files were
made: it keeps the GTF lines of the genes with the most reads in each
compartment, replaces the bases more than 100 bp from these lines with N, and
keeps the alignments that overlap these lines. It does not reproduce the files
exactly. It selects 20
genes in each compartment, not 309 genes, and its BAM files have no alignments
outside the genes. It reads the complete files, which are not in this
repository.

## Attribution and terms of use

These files are not covered by the GPL licence of the RiboseQC code. They stay
under the terms of their sources:

- **Reads.** The data are open access in the NCBI GEO and SRA. NCBI places no
  restriction on their use or distribution. The submitters can claim rights in
  the data, and NCBI cannot grant permission for them. The BAM files hold a
  small subset of the alignments of the study, with attribution to the
  authors [1].
- **Genome sequence and annotation.** TAIR distributes the TAIR10 genome
  sequence and the Araport11 annotation under its terms of use
  (https://phoenixbioinformatics.atlassian.net/wiki/display/COM/TAIR+Terms+of+Use).
  They allow you to create works based on the data and to redistribute these
  works freely. They ask for attribution to TAIR: cite at least one TAIR
  publication [4] and give the URL of the terms of use.

If you use these files in a publication, cite the papers below.

## References

1. Hsu PY, Calviello L, Wu HL, Li FW, Rothfels CJ, Ohler U, Benfey PN (2016).
   Super-resolution ribosome profiling reveals unannotated translation events
   in *Arabidopsis*. *PNAS* 113, E7126–E7135. doi:10.1073/pnas.1614788113
2. Cheng CY, Krishnakumar V, Chan AP, Thibaud-Nissen F, Schobel S, Town CD
   (2017). Araport11: a complete reannotation of the *Arabidopsis thaliana*
   reference genome. *The Plant Journal* 89, 789–804. doi:10.1111/tpj.13415
3. Arabidopsis Genome Initiative (2000). Analysis of the genome sequence of the
   flowering plant *Arabidopsis thaliana*. *Nature* 408, 796–815.
   doi:10.1038/35048692
4. Lamesch P et al. (2012). The Arabidopsis Information Resource (TAIR):
   improved gene annotation and new tools. *Nucleic Acids Research* 40,
   D1202–D1210. doi:10.1093/nar/gkr1090

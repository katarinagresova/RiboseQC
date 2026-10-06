# Test data

A small data set for the tests of `prepare_annotation_files()`. It covers 7
genes of human chr22 (GRCh38): *TXNRD2*, *RANBP1*, *PPIL2*, *TOP3B*, *DDT*,
*TCN2* and *PISD*. The files are copies of ORFquant's example data
(`inst/extdata/` of ORFquant, which has the script that made them,
`inst/scripts/make_example_data.R`).

| File | Content | Source |
|---|---|---|
| `chr22_example.gtf.gz` | All GENCODE release 47 lines of the 7 genes, unchanged. | GENCODE [1] |
| `chr22_example.fa.gz`, `.fai`, `.gzi` | The sequence of chr22 from the GRCh38 primary assembly (GENCODE's `GRCh38.primary_assembly.genome.fa`). Every base more than 1 kb from the 7 genes is replaced by N. Compressed with `bgzip` and indexed. | GRCh38, Genome Reference Consortium [2] |

## Attribution and terms of use

These files are not covered by the GPL licence of the RiboseQC code. They stay
under the terms of their sources. GENCODE data are open access. EMBL-EBI
places no restriction on their use or redistribution other than those of the
original data owners, and expects attribution [1]. None of the sources states
a licence. If you use these files in a publication, cite the papers below.

## References

1. Mudge JM et al. (2025). GENCODE 2025: reference gene annotation for human
   and mouse. *Nucleic Acids Research* 53, D966–D975.
   doi:10.1093/nar/gkae1078
2. Schneider VA et al. (2017). Evaluation of GRCh38 and de novo haploid genome
   assemblies demonstrates the enduring quality of the reference assembly.
   *Genome Research* 27, 849–864. doi:10.1101/gr.213611.116

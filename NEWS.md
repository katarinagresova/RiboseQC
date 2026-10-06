# RiboseQC 1.2.0

## Changes in results

- `RiboseQC_analysis()` now gives the same P-site offsets in each run on the
  same data. Before, the offsets of some read lengths could change from one
  run to the next: `calc_cutoffs_from_profiles()` and `choose_readlengths()`
  call `kmeans()`, which picks random start centres, and RiboseQC set no
  seed. Each `kmeans()` call now uses a fixed seed (666) and R's default kinds
  of random numbers, through the withr package (which ggplot2 already
  installs). The random numbers of your R session don't change. So the
  offsets of these read lengths, and the P-sites, bedgraphs, `*_for_ORFquant`
  files and report plots that use them, can differ from those of an earlier
  run. The offsets of read lengths with a clear frame signal don't change.

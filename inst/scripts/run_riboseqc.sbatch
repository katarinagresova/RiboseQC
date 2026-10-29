#!/bin/bash
## SLURM job running RiboseQC on one or more samples with run_riboseqc.R,
## installed with the package next to this file. Its arguments are passed to
## run_riboseqc.R, with --cores set to the job's CPUs:
##
##   sbatch run_riboseqc.sbatch --gtf genes.gtf --fasta genome.fa \
##     --bam sample.bam --outdir results
##
## Submit it from a shell in which Rscript is the R with RiboseQC (e.g. after
## conda activate riboseqc): the job inherits the shell's environment. Change
## the resources below on the sbatch command line (e.g. sbatch -c 4 --mem=64G)
## and add --partition or --account if your cluster needs them.
##
## Resources, measured on human samples with GENCODE 47, the whole genome, with
## one CPU, the report and rescue_all_rls = TRUE: the annotation took 9 min and
## 7 GB of memory; RiboseQC_analysis() took 4.5 min and 2.5 GB for a sample
## with a 69 MB BAM, and 23 min and 18 GB for one with a 2.3 GB BAM. Time and
## memory grow with the size of the BAM, and the time depends on the CPU. The
## samples of one job run one after the other, so their times add up. The 32 GB
## below are enough for a sample with a 2.3 GB BAM, and the 4 hours for the
## annotation and several such samples: ask for more for larger BAMs or more
## samples. With more CPUs, RiboseQC_analysis() reads more chunks of
## alignments at a time, each in a forked copy of R, so memory grows with the
## CPUs: for the 2.3 GB BAM, 2 CPUs took 18 min and 24 GB, and 4 CPUs 16 min
## and 35 GB (ask for about 48 GB). Only the passes over the BAM use
## more than one CPU, not the profiles, the output files and the report. The
## annotation also writes temporary files, in TMPDIR (by default /tmp): if the
## job's /tmp is too small, the annotation stops with "database or disk is
## full"; ask for more temporary space (--tmp, or what your cluster uses) or
## set TMPDIR to a larger directory.
##
## For more samples, give them together to --bam, or reuse the first job's
## annotation once it is done:
##
##   sbatch --dependency=afterok:<its job id> run_riboseqc.sbatch \
##     --annotation results/annotation/genes.gtf_Rannot \
##     --bam sample2.bam --outdir results
#SBATCH --job-name=riboseqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail
script=$(Rscript -e 'cat(system.file("scripts", "run_riboseqc.R", package = "RiboseQC"))')
if [ -z "$script" ]; then
  echo "Error: RiboseQC is not installed in the R of $(command -v Rscript)" >&2
  exit 1
fi
Rscript "$script" --cores "$SLURM_CPUS_PER_TASK" "$@"

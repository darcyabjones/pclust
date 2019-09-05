#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(DECIPHER)))
data(PFASUM)

fas <- "/dev/stdin"
seqs <- readAAStringSet(fas)

# This reduced alphabet is taken from https://doi.org/10.1002/prot.24936
REDUCED <- c("FWY", "ILMV", "C", "DE", "K", "GNQS", "PT", "A", "HR")
aligned <- AlignSeqs(
  seqs,
  alphabet=REDUCED,
  iterations=2,
  refinements=2,
  substitutionMatrix=PFASUM[,,15],
  verbose=FALSE
)

writeXStringSet(aligned, "/dev/stdout")


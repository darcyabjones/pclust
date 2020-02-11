#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(DECIPHER)))
data(PFASUM)

fas <- "/dev/stdin"
seqs <- readAAStringSet(fas)
seqs <- RemoveGaps(seqs, removeGaps = "all", processors = 1)

# This reduced alphabet is taken from https://doi.org/10.1002/prot.24936
#REDUCED <- c("FWY", "ILMV", "C", "DE", "K", "GNQS", "PT", "A", "HR")
#  alphabet=REDUCED,
aligned <- AlignSeqs(
  seqs,
  iterations=2,
  refinements=2,
  substitutionMatrix=PFASUM[,,15],
  processors=1,
  verbose=FALSE
)

adjusted <- AdjustAlignment(aligned, substitutionMatrix = PFASUM[,,15], processors = 1)
staggered <- StaggerAlignment(adjusted, processors = 1, verbose = FALSE)
consensus <- ConsensusSequence(
  staggered,
  threshold = 0.50,
  ignoreNonBases = FALSE,
  ambiguity = FALSE,
  includeTerminalGaps = TRUE,
  noConsensusChar = "-"
)

names(consensus) <- paste(names(staggered)[1], "_consensus", sep = '')

staggered_cons <- c(consensus, staggered)

writeXStringSet(staggered_cons, "/dev/stdout")

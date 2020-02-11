#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(DECIPHER)))

fas <- "/dev/stdin"
seqs <- readAAStringSet(fas)

consensus <- ConsensusSequence(
  seqs,
  threshold = 0.50,
  ignoreNonBases = FALSE,
  ambiguity = FALSE,
  includeTerminalGaps = TRUE,
  noConsensusChar = "-"
)

names(consensus) <- paste(names(staggered)[1], "_consensus", sep = '')

seqs_cons <- c(consensus, seqs)

writeXStringSet(seqs_cons, "/dev/stdout")

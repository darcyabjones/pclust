#!/usr/bin/env Rscript

VERSION = "0.0.1"

suppressPackageStartupMessages(library("optparse"))
suppressWarnings(suppressPackageStartupMessages(library(DECIPHER)))
data(PFASUM)


option_list <- list(
    make_option(
        c("-i", "--infile"),
        type="character",
        action="store",
        help="The input fasta to align (required)."
    ),
    make_option(
        c("-o", "--outfile"),
        type="character",
        action="store",
        help="The output file to write to (required)."
    ),
    make_option(
        c("-c", "--add-consensus"),
        type="logical",
        dest="consensus",
        action="store_true",
        default=FALSE,
        help="Add a consensus sequence to the MSA."
    ),
    make_option(
        c("-s", "--no-stagger"),
        type="logical",
        dest="stagger",
        action="store_false",
        default=TRUE,
        help="Don't do the staggering step."
    ),
    make_option(
      c("-m", "--matrix"),
      type="integer",
      action="store",
      default=15,
      help=paste(
        "Which PFASUM substitution matrix should we use.",
        "11-100 are available.",
        "Default 15."
      )
    ),
    make_option(
      c("-t", "--consensus-threshold"),
      type="double",
      dest="consensus_threshold",
      action="store",
      default=0.5,
      help=paste(
        "The threshold for proportion of sequences represented in a column for consensus finding.",
        "Must be >0 and < 1.",
        "Default 0.5 ."
      )
    ),
    make_option(
        c("-v", "--verbose"),
        type="logical",
        action="store_true",
        default=FALSE,
        help="Print logging info to stderr.",
    ),
    make_option(
      "--threads",
      type="integer",
      action="store",
      default=1,
      help="How many threads can decipher use? Default: 1"
    ),
    make_option(
        "--version",
        type="logical",
        action="store_true",
        default=FALSE,
        help="Print version and exit.",
    )
)

parser <- OptionParser(
    usage = "%prog --infile in.fasta --outfile out.fasta",
    option_list = option_list
)

args <- parse_args(parser)

log_stderr <- function(...) {
  cat(sprintf(...), sep='', file=stderr())
}

quit_with_err <- function(...) {
  log_stderr(...)
  quit(save = "no", status = 1, runLast = FALSE)
}

validate_file <- function(path) {
  if (is.null(path)) {
    quit_with_err("Please provide required file\n")
  }
}


align_seqs <- function(
  seqs,
  stagger=TRUE,
  sub_matrix=30,
  threads=1,
  verbose=FALSE
) {
  subs_matrix <- PFASUM[,,as.character(sub_matrix)]
  aligned <- AlignSeqs(
    seqs,
    iterations=2,
    refinements=2,
    substitutionMatrix=subs_matrix,
    processors=threads,
    verbose=verbose
  )

  if (verbose) {
    log_stderr("aligned seqs\n")
  }

  adjusted <- AdjustAlignment(
    aligned,
    substitutionMatrix = subs_matrix,
    processors = threads
  )

  if (verbose) {
    log_stderr("adjusted alignment\n")
  }

  if (stagger) {
    staggered <- StaggerAlignment(
      adjusted,
      processors = threads,
      verbose = verbose
    )
    if (verbose) {
      log_stderr("staggered alignment\n")
    }
  } else {
    staggered <- adjusted
  }

  return(staggered)
}


add_consensus <- function(aligned, threshold) {
    consensus <- ConsensusSequence(
      aligned,
      threshold = threshold,
      ignoreNonBases = FALSE,
      ambiguity = FALSE,
      includeTerminalGaps = TRUE,
      noConsensusChar = "-"
    )

    names(consensus) <- paste0(
      gsub("(\\S*).*", "\\1", names(aligned)[1], perl=TRUE),
      "_consensus"
    )
    with_consensus <- c(consensus, aligned)

    return(with_consensus)
}


main <- function(args) {
  if (args$version) {
    cat(VERSION, file=stdout())
    quit(save = "no", status = 0, runLast = FALSE)
  }

  # anything that normally goes to stdout will go to stderr.
  sink(stderr(), type = "output")
  if (args$verbose) {
    log_stderr("start\n")
  }

  validate_file(args$infile)
  validate_file(args$outfile)

  if ((args$matrix < 11) || (args$matrix > 100)) {
    quit_with_err("The matrix value must be between 11 and 100 (inclusive).")
  }

  if ((args$consensus_threshold <= 0) || (args$consensus_threshold >= 1)) {
    quit_with_err("The consensus threshold value must be between 0 and 1 (exclusive).")
  }

  seqs <- readAAStringSet(args$infile)
  seqs <- RemoveGaps(seqs, removeGaps = "all", processors = args$threads)

  log_stderr(as.character(length(seqs)))

  if (args$verbose) {
    log_stderr("read file and removed gaps\n")
  }

  if (length(seqs) > 1) {
    aligned <- align_seqs(
      seqs,
      args$stagger,
      args$matrix,
      args$threads,
      args$verbose
    )
  } else if (length(seqs) == 0) {
    quit_with_err("Input file was empty.\n")
  } else {
    aligned <- seqs
  }

  if (args$consensus) {
    with_consensus <- add_consensus(aligned, args$consensus_threshold)

    if (args$verbose) {
      log_stderr("added consensus sequence\n")
    }

  } else {
    with_consensus <- aligned
  }

  writeXStringSet(with_consensus, args$outfile)

  if (args$verbose) {
    log_stderr("Finished!\n")
  }
}

main(args)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library("readr")
    library("tibble")
    library("dplyr")
    library("magrittr")
    library("stringr")
})


usage <- function() {
    return("add_attribute_col.R <gff3>")
}

args <- commandArgs(trailingOnly=TRUE)

# If no args provided, kill it.
if (length(args) != 1) {
    print(usage())
    quit(save = "no", status = 1, runLast = FALSE)
}

gff_path <- args[1]

read_gff <- function(path) {
    col_names <- c(
        "seqid",
        "analysis",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    )

    col_types <- cols(
        seqid = col_character(),
        analysis = col_character(),
        type = col_character(),
        start = col_integer(),
        end = col_integer(),
        score = col_character(),
        strand = col_character(),
        phase = col_character(),
        attributes = col_character()
    )
    return(read_tsv(path, col_names = col_names, col_type = col_types))
}

gff <- read_gff(gff_path) %>%
    mutate(
        attributes = {
            ids <- str_match(attributes, "ID=([^;]*)")[,2]
            ids[!is.na(ids)] <- paste0("original_id=", ids[!is.na(ids)], ";")
            ids[is.na(ids)] <- ""
            parents <- str_match(attributes, "Parent=([^;]*)")[,2]
            parents[!is.na(parents)] <- paste0("original_parent=",
                                               parents[!is.na(parents)], ";")
            parents[is.na(parents)] <- ""
            paste0(ids, parents, attributes)
        }
    )

# Write to stdout
cat(format_tsv(gff, col_names = FALSE))

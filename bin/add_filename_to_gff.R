#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library("readr")
    library("tibble")
    library("dplyr")
    library("magrittr")
    library("stringr")
})


usage <- function() {
    return("add_filename_to_gff.R <gff3> <out_gff3> <out_map>")
}

args <- commandArgs(trailingOnly=TRUE)

# If no args provided, kill it.
if (length(args) != 3) {
    print(usage())
    quit(save = "no", status = 1, runLast = FALSE)
}

gff_path <- args[1]
gff_out_path <- args[2]
out_map_path <- args[3]

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
    return(read_tsv(path, col_names = col_names, col_type = col_types, comment="#"))
}

base <- basename(tools::file_path_sans_ext(gff_path))

gff <- read_gff(gff_path) %>%
    mutate(
        attributes = gsub("ID=", paste0("ID=", base, "."), attributes) %>%
            gsub("Parent=", paste0("Parent=", base, "."), .)
    ) %>%
    write_tsv(gff_out_path, col_names = FALSE)


gff %>%
    select(type, attributes) %>%
    mutate(
        genome = base,
        id = str_match(attributes, "ID=([^;]*)")[,2],
        original_id = str_match(attributes, "original_id=([^;]*)")[,2],
    ) %>%
    select(genome, type, id, original_id) %>%
    filter(!(is.na(id) && is.na(original_id))) %>%
    write_tsv(out_map_path, col_names = TRUE)

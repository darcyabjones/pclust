#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library("readr")
    library("tibble")
    library("dplyr")
    library("magrittr")
})


usage <- function() {
    return("join_clusters.R <dedup> <cascade> <profile> <profile_stats>")
}

args <- commandArgs(trailingOnly=TRUE)

# If no args provided, kill it.
if (length(args) != 4) {
    print(usage())
    quit(save = "no", status = 1, runLast = FALSE)
}

read_cluster <- function(path) {
    coltypes <- cols(
        cluster = col_character(),
        member = col_character()
    )

    return(read_tsv(path, col_type = coltypes))
}


read_clusters <- function(dedup, cascade, profile) {
    dedup <- read_cluster(dedup) %>%
        rename(dedup = cluster, gene = member)

    cascade <- read_cluster(cascade) %>%
        rename(cascade = cluster, dedup = member)

    profile <- read_cluster(profile) %>%
        rename(profile = cluster, dedup = member)

    joined <- dedup %>%
        full_join(cascade, by="dedup") %>%
        full_join(profile, by="dedup")

    return(joined)
}


count_clusters <- function(df) {
    dedup_count <- df %>%
        group_by(dedup) %>%
        summarize(dedup_count = n())

    cascade_count <- df %>%
        group_by(cascade) %>%
        summarize(cascade_count = n())

    profile_count <- df %>%
        group_by(profile) %>%
        summarize(profile_count = n())

    joined <- df %>%
        left_join(dedup_count, by="dedup") %>%
        left_join(cascade_count, by="cascade") %>%
        left_join(profile_count, by="profile") %>%
        select(gene, dedup, dedup_count, cascade,
               cascade_count, profile, profile_count)

    return(joined)
}


read_cluster_stat <- function(path) {
    coltypes <- cols(
        query = col_character(),
        target = col_character(),
        ident = col_double(),
        length = col_integer(),
        mismatch = col_integer(),
        ngap = col_integer(),
        qstart = col_integer(),
        qend = col_integer(),
        tstart = col_integer(),
        tend = col_integer(),
        evalue = col_double(),
        bitscore = col_integer()
    )
    return(read_tsv(path, col_type = coltypes))
}


read_cluster_stats <- function(clusters, path) {
    stats <- read_cluster_stat(path) %>%
        rename(dedup = target) %>%
        select(-query)


    joined <- clusters %>%
        left_join(stats, by="dedup")

    return(joined)
}

joined <- read_clusters(args[1], args[2], args[3]) %>%
    count_clusters() %>%
    read_cluster_stats(args[4]) %>%
    arrange(desc(profile_count), profile,
            desc(cascade_count), cascade,
            desc(ident),
            desc(dedup_count), dedup,
            gene)

# Write to stdout
cat(format_tsv(joined))

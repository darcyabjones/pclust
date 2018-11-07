#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library("readr")
    library("tibble")
    library("dplyr")
    library("magrittr")
    library("purrr")
})


usage <- function() {
    return("summarise_annotation_clusters.R <clusters> <annotations> <joined> <summarised>")
}

args <- commandArgs(trailingOnly=TRUE)

# If no args provided, kill it.
if (length(args) != 4) {
    print(usage())
    quit(save = "no", status = 1, runLast = FALSE)
}

clusters_path <- args[1]
annotations_path <- args[2]
joined_path <- args[3]
summarized_path <- args[4]

annotations_cols <- cols(
    seqid = col_character(),
    prediction = col_character(),
    score = col_double(),
    analysis = col_character(),
    start = col_integer(),
    end = col_integer(),
    attributes = col_character()
)
annotations <- read_tsv(annotations_path, col_types = annotations_cols)


clusters_cols <- cols(
  gene = col_character(),
  dedup = col_character(),
  dedup_count = col_integer(),
  cascade = col_character(),
  cascade_count = col_integer(),
  profile = col_character(),
  profile_count = col_integer(),
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
clusters <- read_tsv(clusters_path, col_types = clusters_cols)

joined <- full_join(
    select(clusters, gene, dedup, cascade, profile),
    annotations,
    by = c("dedup"="seqid")
)

write_tsv(joined, joined_path)


# Summarise clusters

size <- clusters %>%
    group_by(profile) %>%
    summarize(cluster_size = n(), cluster_length = mean(length, na.rm = TRUE))


effectorp <- joined %>%
    filter(analysis == "effectorp") %>%
    group_by(profile) %>%
    summarize(
        n_effector = sum(prediction == "effector", na.rm = TRUE),
        p_effector = mean(prediction == "effector", na.rm = TRUE),
        n_unlikely_effector = sum(prediction == "unlikely effector", na.rm = TRUE),
        p_unlikely_effector = mean(prediction == "unlikely effector", na.rm = TRUE),
    ) %>%
    rename_at(vars(-profile), ~ paste0("effectorp_", .))


apoplastp <- joined %>%
    filter(analysis == "apoplastp") %>%
    group_by(profile) %>%
    summarize(
        n_apoplastic = sum(prediction == "apoplastic", na.rm = TRUE),
        p_apoplastic = mean(prediction == "apoplastic", na.rm = TRUE)
    ) %>%
    rename_at(vars(-profile), ~ paste0("apoplastp_", .))

localizer_plant <- joined %>%
    filter(analysis == "localizer_plant") %>%
    group_by(profile) %>%
    summarize(
        n_nucleus = sum(prediction == "nucleus", na.rm = TRUE),
        p_nucleus = mean(prediction == "nucleus", na.rm = TRUE),
        n_mitochondria = sum(prediction == "mitochondria", na.rm = TRUE),
        p_mitochondria = mean(prediction == "mitochondria", na.rm = TRUE),
        n_chloroplast = sum(prediction == "chloroplast", na.rm = TRUE),
        p_chloroplast = mean(prediction == "chloroplast", na.rm = TRUE),
    ) %>%
    rename_at(vars(-profile), ~ paste0("localizer_plant_", .))


# TODO: Modify location to be original location based on signalp4 output.
localizer_effector <- joined %>%
    filter(analysis == "localizer_effector") %>%
    group_by(profile) %>%
    summarize(
        n_nucleus = sum(prediction == "nucleus", na.rm = TRUE),
        p_nucleus = mean(prediction == "nucleus", na.rm = TRUE),
        n_mitochondria = sum(prediction == "mitochondria", na.rm = TRUE),
        p_mitochondria = mean(prediction == "mitochondria", na.rm = TRUE),
        n_chloroplast = sum(prediction == "chloroplast", na.rm = TRUE),
        p_chloroplast = mean(prediction == "chloroplast", na.rm = TRUE),
    ) %>%
    rename_at(vars(-profile), ~ paste0("localizer_effector_", .))

targetp_non_plant <- joined %>%
    filter(analysis == "targetp_non_plant", score > 0) %>%
    group_by(profile) %>%
    summarize(
        n_signal = sum(prediction == "secreted", na.rm = TRUE),
        p_signal = mean(prediction == "secreted", na.rm = TRUE),
        n_mitochondria = sum(prediction == "mitochondria", na.rm = TRUE),
        p_mitochondria = mean(prediction == "mitochondria", na.rm = TRUE),
    ) %>%
    rename_at(vars(-profile), ~ paste0("targetp_non_plant_", .))

targetp_plant <- joined %>%
    filter(analysis == "targetp_plant", score > 0) %>%
    group_by(profile) %>%
    summarize(
        n_signal = sum(prediction == "secreted", na.rm = TRUE),
        p_signal = mean(prediction == "secreted", na.rm = TRUE),
        n_mitochondria = sum(prediction == "mitochondria", na.rm = TRUE),
        p_mitochondria = mean(prediction == "mitochondria", na.rm = TRUE),
        n_chloroplast = sum(prediction == "chloroplast", na.rm = TRUE),
        p_chloroplast = mean(prediction == "chloroplast", na.rm = TRUE),
    ) %>%
    rename_at(vars(-profile), ~ paste0("targetp_plant_", .))

signalp3_hmm <- joined %>%
    filter(analysis == "signalp3_hmm") %>%
    group_by(profile) %>%
    summarize(
        n_signal = sum(prediction == "signal_peptide", na.rm = TRUE),
        p_signal = mean(prediction == "signal_peptide", na.rm = TRUE),
    ) %>%
    rename_at(vars(-profile), ~ paste0("signalp3_hmm_", .))

signalp4 <- joined %>%
    filter(analysis == "signalp4") %>%
    group_by(profile) %>%
    summarize(
        n_signal = sum(prediction == "signal_peptide", na.rm = TRUE),
        p_signal = mean(prediction == "signal_peptide", na.rm = TRUE),
    ) %>%
    rename_at(vars(-profile), ~ paste0("signalp4_", .))


phobius_sec <- joined %>%
    filter(analysis == "phobius") %>%
    group_by(profile) %>%
    summarize(
        n_signal = sum(prediction == "signal_peptide", na.rm = TRUE),
        p_signal = mean(prediction == "signal_peptide", na.rm = TRUE),
    ) %>%
    rename_at(vars(-profile), ~ paste0("phobius_", .))


phobius_tm <- joined %>%
    filter(analysis == "phobius") %>%
    group_by(gene) %>%
    summarize(
        profile = first(profile),
        n_tm = sum(prediction == "transmembrane", na.rm = TRUE),
        len_tm = sum((end - start + 1)[prediction == "transmembrane"], na.rm = TRUE),
        n_apo = sum(prediction == "non_cytoplasmic", na.rm = TRUE),
        len_apo = sum((end - start + 1)[prediction == "non_cytoplasmic"], na.rm = TRUE),
        n_cyto = sum(prediction == "cytoplasmic", na.rm = TRUE),
        len_cyto = sum((end - start + 1)[prediction == "cytoplasmic"], na.rm = TRUE),
    ) %>%
    group_by(profile) %>%
    summarize(
        n_transmembrane = mean(n_tm, na.rm = TRUE),
        p_transmembrane = mean(n_tm > 0, na.rm = TRUE),
        l_transmembrane = mean(len_tm, na.rm = TRUE),
        n_cytoplasmic = mean(n_cyto, na.rm = TRUE),
        p_cytoplasmic = mean(n_cyto > 0, na.rm = TRUE),
        l_cytoplasmic = mean(len_cyto, na.rm = TRUE),
        n_non_cytoplasmic = mean(n_apo, na.rm = TRUE),
        p_non_cytoplasmic = mean(n_apo > 0, na.rm = TRUE),
        l_non_cytoplasmic = mean(len_apo, na.rm = TRUE),
    ) %>%
    rename_at(vars(-profile), ~ paste0("phobius_", .))


tmhmm <- joined %>%
    filter(analysis == "tmhmm") %>%
    group_by(gene) %>%
    summarize(
        profile = first(profile),
        n_tm = sum(prediction == "transmembrane", na.rm = TRUE),
        len_tm = sum((end - start + 1)[prediction == "transmembrane"], na.rm = TRUE),
        n_apo = sum(prediction == "non_cytoplasmic", na.rm = TRUE),
        len_apo = sum((end - start + 1)[prediction == "non_cytoplasmic"], na.rm = TRUE),
        n_cyto = sum(prediction == "cytoplasmic", na.rm = TRUE),
        len_cyto = sum((end - start + 1)[prediction == "cytoplasmic"], na.rm = TRUE),
    ) %>%
    group_by(profile) %>%
    summarize(
        n_transmembrane = mean(n_tm, na.rm = TRUE),
        p_transmembrane = mean(n_tm > 0, na.rm = TRUE),
        l_transmembrane = mean(len_tm, na.rm = TRUE),
        n_cytoplasmic = mean(n_cyto, na.rm = TRUE),
        p_cytoplasmic = mean(n_cyto > 0, na.rm = TRUE),
        l_cytoplasmic = mean(len_cyto, na.rm = TRUE),
        n_non_cytoplasmic = mean(n_apo, na.rm = TRUE),
        p_non_cytoplasmic = mean(n_apo > 0, na.rm = TRUE),
        l_non_cytoplasmic = mean(len_apo, na.rm = TRUE),
    ) %>%
    rename_at(vars(-profile), ~ paste0("tmhmm_", .))


output <- purrr::reduce(
    list(
        size,
        apoplastp,
        effectorp,
        signalp3_hmm,
        signalp4,
        phobius_sec,
        phobius_tm,
        tmhmm,
        targetp_non_plant,
        targetp_plant,
        localizer_effector,
        localizer_plant
    ),
    function(x, y) full_join(x, y, by="profile")
)

write_tsv(output, summarized_path)

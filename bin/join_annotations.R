#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library("readr")
    library("tibble")
    library("dplyr")
    library("magrittr")
    library("stringr")
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

process_loc <- function(df) {
    regex <- "Y \\((\\d\\.\\d+)\\s+\\|\\s+(\\d+)-(\\d+)\\)"

    chloro <- df %>% select(seqid, chloroplast) %>% filter(chloroplast != "-")
    chloro_matches <- str_match(chloro$chloroplast, regex)
    colnames(chloro_matches) <- c("junk", "score", "start", "end")
    rownames(chloro_matches) <- chloro$seqid
    chloro <- as.tibble(chloro_matches, rownames = "seqid") %>%
        type_convert(
            col_types = cols(
                seqid = col_character(),
                junk = col_character(),
                score = col_double(),
                start = col_integer(),
                end = col_integer()
            )
        ) %>%
        select(-junk) %>%
        mutate(prediction = "chloroplast")

    mito <- df %>% select(seqid, mitochondria) %>% filter(mitochondria != "-")
    mito_matches <- str_match(mito$mitochondria, regex)
    colnames(mito_matches) <- c("junk", "score", "start", "end")
    rownames(mito_matches) <- mito$seqid
    mito <- as.tibble(mito_matches, rownames = "seqid") %>%
        type_convert(
            col_types = cols(
                seqid = col_character(),
                junk = col_character(),
                score = col_double(),
                start = col_integer(),
                end = col_integer()
            )
        ) %>%
        select(-junk) %>%
        mutate(prediction = "mitochondria")

    nucl <- df %>%
        select(seqid, nucleus) %>%
        filter(nucleus != "-") %>%
        mutate(peptide = str_match(nucleus, "\\((.+)\\)")[,2]) %>%
        select(-nucleus) %>%
        separate_rows(peptide, sep = ",") %>%
        mutate(prediction = "nucleus", attributes = paste0("peptide=", peptide)) %>%
        select(-peptide)

    combined <- bind_rows(chloro, mito, nucl) %>% mutate(analysis = "LOCALIZER")
    return(combined)
}


process_phobius <- function(df) {
    signal <- df %>%
        filter(sp == "Y") %>%
        mutate(cleavage = str_match(prediction, "(\\d+)/\\d+")[, 2],
               nstart = str_match(prediction, "n(\\d+)-")[, 2],
               cstart = str_match(prediction, "-(\\d+)c")[, 2]
               )

    sp <- signal %>%
        mutate(start = as.integer(1),
               end = as.integer(cleavage),
               prediction = "signal_peptide") %>%
        select(seqid, start, end, prediction)

    n_region <- signal %>%
        mutate(start = as.integer(1),
               end = as.integer(nstart) - as.integer(1),
               prediction = "signal_peptide_n_region") %>%
        select(seqid, start, end, prediction)

    h_region <- signal %>%
        mutate(start = as.integer(nstart),
               end = as.integer(cstart),
               prediction = "signal_peptide_h_region") %>%
        select(seqid, start, end, prediction)

    c_region <- signal %>%
        mutate(start = as.integer(cstart) + as.integer(1),
               end = as.integer(cleavage),
               prediction = "signal_peptide_c_region") %>%
        select(seqid, start, end, prediction)

    tm_domains <- df %>%
        filter(tm > 0) %>%
        select(seqid, prediction)

    doms <- str_match_all(tm_domains$prediction, "[io](\\d+)-(\\d+)") %>%
        lapply(function(x) {
            tibble(start = as.integer(x[, 2]),
                   end = as.integer(x[, 3]))
            }
        ) %>%
        mapply(
        function(x, y) {x["seqid"] = y; x},
        .,
        tm_domains$seqid,
        SIMPLIFY=FALSE
        ) %>%
        bind_rows() %>%
        mutate(prediction = "transmembrane")

    combined <- bind_rows(sp, n_region, h_region, c_region, doms) %>%
        mutate(analysis = "phobius")
    return(combined)
}


process_tmhmm <- function(df) {
    tm_domains <- str_match_all(df$topology, "[io](\\d+)-(\\d+)") %>%
        lapply(function(x) {
            tibble(start = as.integer(x[, 2]),
                   end = as.integer(x[, 3]))
            }
        ) %>%
        mapply(
        function(x, y) {x["seqid"] = y; x},
        .,
        df$seqid,
        SIMPLIFY=FALSE
        ) %>%
        bind_rows() %>%
        mutate(prediction = "transmembrane")

    cyto_domains <- str_match_all(df$topology, "(\\d+)i(\\d+)") %>%
        lapply(function(x) {
            tibble(start = as.integer(x[, 2]) + as.integer(1),
                   end = as.integer(x[, 3]) - as.integer(1))
            }
        ) %>%
        mapply(
            function(x, y) {
                if (nrow(x) > 0) {
                    x["seqid"] = y
                    return(x)
                } else {
                    return(x)
                }
            },
            .,
            df$seqid,
            SIMPLIFY=FALSE
        ) %>%
        bind_rows() %>%
        mutate(prediction = "cytoplasmic")

    apo_domains <- str_match_all(df$topology, "(\\d+)o(\\d+)") %>%
        lapply(function(x) {
            tibble(start = as.integer(x[, 2]) + as.integer(1),
                   end = as.integer(x[, 3]) - as.integer(1))
            }
        ) %>%
        mapply(
            function(x, y) {
                if (nrow(x) > 0) {
                    x["seqid"] = y
                    return(x)
                } else {
                    return(x)
                }
            },
            .,
            df$seqid,
            SIMPLIFY=FALSE
        ) %>%
        bind_rows() %>%
        mutate(prediction = "non_cytoplasmic")

    first_domains <- str_match_all(df$topology, "^[io](\\d+)") %>%
        lapply(function(x) {
            tibble(start = as.integer(1),
                   end = as.integer(x[, 2]) - as.integer(1))
            }
        ) %>%
        mapply(
            function(x, y, z) {
                if (nrow(x) == 0) {
                    return(x)
                }

                x["seqid"] = y

                if (startsWith(z, "i")) {
                    x["prediction"] = "cytoplasmic"
                } else {
                    x["prediction"] = "non_cytoplasmic"
                }
                return(x)
            },
            .,
            df$seqid,
            df$topology,
            SIMPLIFY=FALSE
        ) %>%
        bind_rows()

    last_domains <- str_match_all(df$topology, "(\\d+)[io]$") %>%
        lapply(function(x) {
            tibble(start = as.integer(x[, 2]) + as.integer(1))
            }
        ) %>%
        mapply(
            function(x, y, z, w) {
                if (nrow(x) == 0) {
                    return(x)
                }

                x["seqid"] = y
                x["end"] = w
                if (startsWith(z, "i")) {
                    x["prediction"] = "cytoplasmic"
                } else {
                    x["prediction"] = "non_cytoplasmic"
                }
                return(x)
            },
            .,
            df$seqid,
            df$topology,
            df$len,
            SIMPLIFY=FALSE
        ) %>%
        bind_rows()

    combined <- bind_rows(tm_domains, cyto_domains, apo_domains,
                          first_domains, last_domains) %>%
        mutate(analysis = "tmhmm")
    combined <- left_join(combined, select(df, seqid, attributes), by="seqid")
    return(combined)
}


apoplastp <- readr::read_tsv(
    "annotations/apoplastp.tsv",
    col_types=cols(
        seqid = col_character(),
        prediction = col_character(),
        probability = col_double()
    )
) %>%
    rename(score = probability) %>%
    mutate(analysis = "apoplastp")

effectorp <- readr::read_tsv(
    "annotations/effectorp.tsv",
    col_types = cols(
        seqid = col_character(),
        effector = col_character(),
        probability = col_double()
    )
) %>%
    rename(prediction = effector, score = probability) %>%
    mutate(analysis = "effectorp")

localiser_eff <- readr::read_tsv(
    "annotations/localizer_effector.tsv",
    col_types = cols(
        seqid = col_character(),
        chloroplast = col_character(),
        mitochondria = col_character(),
        nucleus = col_character()
    )
) %>% process_loc()


localiser_plant <- readr::read_tsv(
    "annotations/localizer_plant.tsv",
    col_types = cols(
        seqid = col_character(),
        chloroplast = col_character(),
        mitochondria = col_character(),
        nucleus = col_character()
    )
) %>% process_loc()

phobius <- readr::read_tsv(
    "annotations/phobius.tsv",
    col_types = cols(
        seqid = col_character(),
        tm = col_integer(),
        sp = col_character(),
        prediction = col_character()
    )
) %>% process_phobius()

signalp3_hmm <- readr::read_tsv(
    "annotations/signalp3_hmm.tsv",
    col_type = cols(
        seqid = col_character(),
        secreted = col_character(),
        cmax = col_double(),
        pos = col_integer(),
        pos_decision = col_character(),
        sprob = col_double(),
        sprob_decision = col_character()
    )
) %>%
    mutate(start = as.integer(1),
           end = pos,
           score = sprob,
           cmax = paste0("cmax=", cmax),
           sprob = paste0("sprob=", sprob),
           pos = paste0("cmax_pos=", pos),
           attributes = paste(cmax, sprob, pos, sep = ";")
)
signalp3_hmm_conf <- signalp3_hmm %>%
    filter(pos_decision == "Y" & secreted == "S") %>%
    select(seqid, start, end, score, attributes) %>%
    mutate(prediction = "signal_peptide")

signalp3_hmm_low <- signalp3_hmm %>%
    filter(pos_decision != "Y" & secreted == "S") %>%
    select(seqid, score, attributes) %>%
    mutate(prediction = "signal_peptide")

signalp3_hmm_non <- signalp3_hmm %>%
    filter(secreted != "S") %>%
    select(seqid, score, attributes) %>%
    mutate(prediction = "no_signal_peptide")

signalp3_hmm <- bind_rows(signalp3_hmm_conf, signalp3_hmm_low, signalp3_hmm_non) %>%
    mutate(analysis = "signalp3_hmm")

signalp4 <- readr::read_tsv(
    "annotations/signalp4.tsv",
    col_type = cols(
        seqid = col_character(),
        cmax = col_double(),
        cmax_pos = col_integer(),
        ymax = col_double(),
        ymax_pos = col_integer(),
        smax = col_double(),
        smax_pos = col_integer(),
        smean = col_double(),
        d = col_double(),
        secreted = col_character(),
        dmaxcut = col_double(),
        networks = col_character()
    )
) %>%
    mutate(
        score = d,
        start = as.integer(1),
        end = as.integer(ymax_pos - 1),
        cmax = paste0("cmax=", cmax),
        cmax_pos = paste0("cmax_pos=", cmax_pos),
        ymax = paste0("ymax=", ymax),
        ymax_pos = paste0("ymax_pos=", ymax_pos),
        smax = paste0("smax=", smax),
        smax_pos = paste0("smax_pos=", smax_pos),
        smean = paste0("smean=", smean),
        d = paste0("d=", d),
        networks = paste0("networks=", networks),
        attributes = paste(cmax, cmax_pos, ymax, ymax_pos, smax,
                           smax_pos, smean, d, networks, sep = ";")
        )

signalp4_sec <- signalp4 %>%
    filter(secreted == "Y") %>%
    select(seqid, start, end, score, attributes) %>%
    mutate(prediction = "signal_peptide")

signalp4_non <- signalp4 %>%
    filter(secreted == "N") %>%
    select(seqid, score, attributes) %>%
    mutate(prediction = "no_signal_peptide")

signalp4 <- bind_rows(signalp4_sec, signalp4_non) %>%
    mutate(analysis = "signalp4")


targetp_np <- readr::read_tsv(
    "annotations/targetp.tsv",
    col_types = cols(
        seqid = col_character(),
        len = col_integer(),
        mtp = col_double(),
        sp = col_double(),
        other = col_double(),
        loc = col_character(),
        rc = col_integer(),
        tplen = col_character()
    )
) %>% mutate(mtp = paste0("mtp=", mtp),
             sp = paste0("sp=", sp),
             len = paste0("len=", len),
             other = paste0("other=", other),
             score = 1 - ((rc - 1) / 4),
             rc = paste0("rc=", other),
             attributes = paste(mtp, sp, len, other, rc, sep = ";"),
             analysis = "targetp_non_plant"
             )

targetp_pl <- readr::read_tsv(
    "annotations/targetp_plant.tsv",
    col_types = cols(
        seqid = col_character(),
        len = col_integer(),
        ctp = col_double(),
        mtp = col_double(),
        sp = col_double(),
        other = col_double(),
        loc = col_character(),
        rc = col_integer(),
        tplen = col_character()
    )
) %>% mutate(mtp = paste0("mtp=", mtp),
             sp = paste0("sp=", sp),
             ctp = paste0("ctp=", ctp),
             len = paste0("len=", len),
             other = paste0("other=", other),
             score = 1 - ((rc - 1) / 4),
             rc = paste0("rc=", other),
             attributes = paste(mtp, sp, ctp, len, other, rc, sep = ";"),
             analysis = "targetp_plant"
             )

targetp <- bind_rows(targetp_pl, targetp_np) %>%
    mutate(
        prediction = case_when(
            loc == "M" ~ "mitochondria",
            loc == "C" ~ "chloroplast",
            loc == "S" ~ "secreted",
            loc == "_" ~ "other_localisation",
            TRUE ~ "unknown_localisation"
        ),
        start = ifelse(tplen != "-", as.integer(1), NA),
        end = {
            tplen[tplen == "-"] <- NA
            as.integer(tplen)
        },
    ) %>%
    select(seqid, start, end, score, attributes, analysis, prediction)


tmhmm <- readr::read_tsv(
    "annotations/tmhmm.tsv",
    col_types = cols(
        seqid = col_character(),
        len = col_integer(),
        expaa = col_double(),
        first60 = col_double(),
        predhel = col_integer(),
        topology = col_character()
    )
) %>%
    filter(predhel > 0) %>%
    mutate(score = expaa,
           len2 = paste0("len=", len),
           expaa = paste0("expaa=", expaa),
           predhel = paste0("predhel=", predhel),
           first60 = paste0("first60=", first60),
           attributes = paste(len2, expaa, predhel, first60,
                              paste0("topology=", topology), sep = ";")) %>%
    select(seqid, score, len, topology, attributes) %>%
    process_tmhmm()


combined <- bind_rows(effectorp, apoplastp, localiser_eff, localiser_plant, signalp3_hmm, signalp4, targetp, phobius, tmhmm)

# Write to stdout
cat(format_tsv(joined))

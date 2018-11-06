#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library("readr")
    library("tibble")
    library("dplyr")
    library("magrittr")
    library("stringr")
    library("tidyr")
})



#' Extract score from localizer mito and chloro columns.
#' Columns are in the format, (0.999 | 1-20) i.e. prob and start-stop.
#' df must be filtered to contain only positive matches for column.
localizer_scores <- function(df, column) {
    regex <- "Y \\((\\d\\.\\d+)\\s+\\|\\s+(\\d+)-(\\d+)\\)"
    matches <- str_match(pull(df, column), regex)
    colnames(matches) <- c("junk", "score", "start", "end")
    rownames(matches) <- pull(df, "seqid")

    matches <- as.tibble(matches, rownames = "seqid") %>%
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
        mutate(prediction = column)

    return(matches)
}


#' Extract data from localizer output.
process_loc <- function(df) {

    chloro <- df %>%
        select(seqid, chloroplast) %>%
        filter(chloroplast != "-") %>%
        localizer_scores("chloroplast")

    mito <- df %>%
        select(seqid, mitochondria) %>%
        filter(mitochondria != "-") %>%
        localizer_scores("mitochondria")

    nucl <- df %>%
        select(seqid, nucleus) %>%
        filter(nucleus != "-") %>%
        mutate(peptide = str_match(nucleus, "\\((.+)\\)")[,2]) %>%
        select(-nucleus) %>%
        separate_rows(peptide, sep = ",") %>%
        mutate(prediction = "nucleus",
               attributes = paste0("peptide=", peptide)) %>%
        select(-peptide)

    bind_rows(chloro, mito, nucl) %>%
        return()
}


#' Process signalp3 hmm output.
#' This function could be simplified a bit using `case` syntax.
process_signalp3_hmm <- function(df) {
    df <- df %>%
        mutate(start = as.integer(1),
               end = pos,
               score = sprob,
               cmax = paste0("cmax=", cmax),
               sprob = paste0("sprob=", sprob),
               pos = paste0("cmax_pos=", pos),
               attributes = paste(cmax, sprob, pos, sep = ";"))

    # If the cleavage site is confident, include the start, stop.
    conf <- df %>%
        filter(pos_decision == "Y" & secreted == "S") %>%
        select(seqid, start, end, score, attributes) %>%
        mutate(prediction = "signal_peptide")

    # If the cleavage site not confident, dont include start, stop.
    low <- df %>%
        filter(pos_decision != "Y" & secreted == "S") %>%
        select(seqid, score, attributes) %>%
        mutate(prediction = "signal_peptide")

    # If not signal peptide, write different prediction
    non <- df %>%
        filter(secreted != "S") %>%
        select(seqid, score, attributes) %>%
        mutate(prediction = "no_signal_peptide")

    # Combine split dfs and return
    bind_rows(conf, low, non) %>%
        mutate(analysis = "signalp3_hmm") %>%
        return()
}


#' Process signalp4 output.
process_signalp4 <- function(df) {
    df <- df %>%
        mutate(score = d,
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

    sec <- df %>%
        filter(secreted == "Y") %>%
        select(seqid, start, end, score, attributes) %>%
        mutate(prediction = "signal_peptide")

    non <- df %>%
        filter(secreted == "N") %>%
        select(seqid, score, attributes) %>%
        mutate(prediction = "no_signal_peptide")

    bind_rows(sec, non) %>%
        mutate(analysis = "signalp4") %>%
        return()
}


#' Process targetp output, both networks
process_targetp <- function(np, pl) {
    np <- np %>%
        mutate(mtp = paste0("mtp=", mtp),
               sp = paste0("sp=", sp),
               len = paste0("len=", len),
               other = paste0("other=", other),
               score = 1 - ((rc - 1) / 4),
               rc = paste0("rc=", rc),
               attributes = paste(mtp, sp, len, other, rc, sep = ";"),
               analysis = "targetp_non_plant"
               )

    pl <- pl %>%
        mutate(mtp = paste0("mtp=", mtp),
               sp = paste0("sp=", sp),
               ctp = paste0("ctp=", ctp),
               len = paste0("len=", len),
               other = paste0("other=", other),
               score = 1 - ((rc - 1) / 4),
               rc = paste0("rc=", rc),
               attributes = paste(mtp, sp, ctp, len, other, rc, sep = ";"),
               analysis = "targetp_plant"
               )

    bind_rows(pl, np) %>%
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
        select(seqid, start, end, score, attributes, analysis, prediction) %>%
        return()
}


#' Add a seqid column to a dataframe if the the df is not empty.
add_seqid <- function(df, seqid) {
    if (nrow(df) > 0) {
        df["seqid"] = seqid
        return(df)
    } else {
        return(df)
    }
}


#' Generic function to abstract extacting topology info.
#' NB here i'm using mapply as a zip function.
#' FUN should take a matrix as returned from str_match.
extract_domain <- function(topo, seqid, regex, FUN) {
    str_match_all(topo, regex) %>%
        lapply(FUN) %>%
        mapply(add_seqid, ., seqid, SIMPLIFY=FALSE) %>%
        bind_rows() %>%
        return()
}


#' Extract transmembrane domains from standard topology string.
extract_tm_domains <- function(topo, seqid) {
    extract_domain(topo, seqid, "[io](\\d+)-(\\d+)", function(x) {
            tibble(start = as.integer(x[, 2]),
                   end = as.integer(x[, 3]),
                   prediction = "transmembrane")
    })
}


#' Extract regions between transmembrane domains as loops.
extract_loop_domains <- function(topo, seqid) {
    extract_domain(topo, seqid, "(\\d+)([io])(\\d+)", function(x) {
            pred <- ifelse(x[, 3] == 'i', "cytoplasmic", "non_cytoplasmic")
            tibble(start = as.integer(x[, 2]) + as.integer(1),
                   end = as.integer(x[, 4]) - as.integer(1),
                   prediction = pred)
    })
}


#' Extract the leading domain.
extract_first_domain <- function(topo, seqid, start = 1) {
    extract_domain(topo, seqid, "^([io])(\\d+)", function(x) {
        pred <- ifelse(
            startsWith(x[, 2], "i"),
            "cytoplasmic",
            "non_cytoplasmic"
        )

        tibble(start = as.integer(start),
               end = as.integer(x[, 3]) - as.integer(1),
               prediction = pred) %>%
            filter(start != end) %>%
            return()
    })
}


#' Extract the trailing domain.
extract_last_domain <- function(topo, seqid, end) {
    doms <- extract_domain(topo, seqid, "(\\d+)([io])$", function(x) {
        pred <- ifelse(
            startsWith(x[, 3], "i"),
            "cytoplasmic",
            "non_cytoplasmic"
        )

        tibble(start = as.integer(x[, 2]) + as.integer(1), prediction = pred)
    })

    lens <- tibble(seqid = seqid, end = as.integer(end))
    left_join(doms, lens, by = "seqid")
}


process_phobius <- function(df, lengths) {
    df <- df %>%
        mutate(
            cleavage = str_match(prediction, "(\\d+)/\\d")[, 2] %>%
                as.integer() %>% replace_na(as.integer(0)),
            nstart = str_match(prediction, "n(\\d+)-")[, 2],
            cstart = str_match(prediction, "-(\\d+)c")[, 2]
        ) %>%
        left_join(lengths, by = "seqid")

    signal <- df %>%
        filter(sp == "Y") %>%
        mutate(
            start = as.integer(1),
            end = as.integer(cleavage),
            prediction = "signal_peptide"
        ) %>%
        select(seqid, start, end, prediction)

    n_region <- df %>%
        filter(sp == "Y") %>%
        mutate(start = as.integer(1),
               end = as.integer(nstart) - as.integer(1),
               prediction = "signal_peptide_n_region") %>%
        select(seqid, start, end, prediction)

    h_region <- df %>%
        filter(sp == "Y") %>%
        mutate(start = as.integer(nstart),
               end = as.integer(cstart),
               prediction = "signal_peptide_h_region") %>%
        select(seqid, start, end, prediction)

    c_region <- df %>%
        filter(sp == "Y") %>%
        mutate(start = as.integer(cstart) + as.integer(1),
               end = as.integer(cleavage),
               prediction = "signal_peptide_c_region") %>%
        select(seqid, start, end, prediction)

    tm <- filter(df, tm > 0)
    tm_domains <- extract_tm_domains(tm$prediction, tm$seqid)
    loop_domains <- extract_loop_domains(tm$prediction, tm$seqid)
    last_domains <- extract_last_domain(tm$prediction, tm$seqid, tm$len)

    # First domains of secreted proteins are picked up by loop_domains fn.
    non_secreted <- filter(tm, sp != "Y")
    first_domains <- extract_first_domain(
        non_secreted$prediction,
        non_secreted$seqid,
        start = as.integer(1)
    )

    bind_rows(
            signal,
            n_region,
            h_region,
            c_region,
            tm_domains,
            loop_domains,
            first_domains,
            last_domains
        ) %>%
        mutate(analysis = "phobius") %>%
        return()
}



process_tmhmm <- function(df) {
    df <- df %>%
        filter(predhel > 0) %>%
        mutate(
            score = expaa,
            len2 = paste0("len=", len),
            expaa = paste0("expaa=", expaa),
            predhel = paste0("predhel=", predhel),
            first60 = paste0("first60=", first60),
            attributes = paste(
                len2,
                expaa,
                predhel,
                first60,
                paste0("topology=", topology),
                sep = ";"
            )
        ) %>%
        select(seqid, score, len, topology, attributes)

    tm_domains <- extract_tm_domains(df$topology, df$seqid)
    loop_domains <- extract_loop_domains(df$topology, df$seqid)
    first_domains <- extract_first_domain(df$topology, df$seqid)
    last_domains <- extract_last_domain(df$topology, df$seqid, df$len)

    combined <- bind_rows(
            tm_domains,
            loop_domains,
            first_domains,
            last_domains
        ) %>%
        mutate(analysis = "tmhmm") %>%
        left_join(select(df, seqid, attributes), by="seqid") %>%
        return()
}

## Main.

usage <- function() {
    return("join_annotations.R <apoplastp> <effectorp> <localizer_eff> <localizer_plant> <signalp3_hmm> <signalp3_nn> <signalp4> <targtp_non_plant> <targetp_plant> <tmhmm> <phobius>")
}

args <- commandArgs(trailingOnly=TRUE)

# If no args provided, kill it.
if (length(args) != 11) {
    print(usage())
    quit(save = "no", status = 1, runLast = FALSE)
}

apoplastp_path <- args[1]
effectorp_path <- args[2]
localizer_eff_path <- args[3]
localizer_plant_path <- args[4]
signalp3_hmm_path <- args[5]
signalp3_nn_path <- args[6]
signalp4_path <- args[7]
targetp_np_path <- args[8]
targetp_pl_path <- args[9]
tmhmm_path <- args[10]
phobius_path <- args[11]


# Apoplastp
apoplastp_cols <- cols(
    seqid = col_character(),
    prediction = col_character(),
    probability = col_double()
)

apoplastp <- read_tsv(apoplastp_path, col_types = apoplastp_cols) %>%
    rename(score = probability) %>%
    mutate(
        analysis = "apoplastp",
        prediction = str_to_lower(prediction) %>% str_replace_all("-", "_")
    )


# Effectorp
effectorp_cols <- cols(
    seqid = col_character(),
    effector = col_character(),
    probability = col_double()
)

effectorp <- read_tsv(effectorp_path, col_types = effectorp_cols) %>%
    rename(prediction = effector, score = probability) %>%
    mutate(
        analysis = "effectorp",
        prediction = str_to_lower(prediction) %>% str_replace_all("-", "_")
    )


# Localizer
localizer_cols <- cols(
    seqid = col_character(),
    chloroplast = col_character(),
    mitochondria = col_character(),
    nucleus = col_character()
)

localizer_eff <- read_tsv(localizer_eff_path,
                          col_types = localizer_cols) %>%
    process_loc() %>%
    mutate(analysis = "localizer_effector")

localizer_plant <- read_tsv(localizer_plant_path,
                            col_types = localizer_cols) %>%
    process_loc() %>%
    mutate(analysis = "localizer_plant")


# Signalp3
signalp3_hmm_cols <- cols(
    seqid = col_character(),
    secreted = col_character(),
    cmax = col_double(),
    pos = col_integer(),
    pos_decision = col_character(),
    sprob = col_double(),
    sprob_decision = col_character()
)

signalp3_hmm <- read_tsv(signalp3_hmm_path, col_type = signalp3_hmm_cols) %>%
    process_signalp3_hmm()


# Signalp4
signalp4_cols <- cols(
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

signalp4 <- read_tsv(signalp4_path, col_type = signalp4_cols) %>%
    process_signalp4()


# Targetp
targetp_np_cols <- cols(
    seqid = col_character(),
    len = col_integer(),
    mtp = col_double(),
    sp = col_double(),
    other = col_double(),
    loc = col_character(),
    rc = col_integer(),
    tplen = col_character()
)

targetp_pl_cols <- cols(
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


targetp <- process_targetp(
    read_tsv(targetp_np_path, col_types = targetp_np_cols),
    read_tsv(targetp_pl_path, col_types = targetp_pl_cols)
)


# TMHMM
tmhmm_cols <- cols(
    seqid = col_character(),
    len = col_integer(),
    expaa = col_double(),
    first60 = col_double(),
    predhel = col_integer(),
    topology = col_character()
)

tmhmm_tmp <- readr::read_tsv(tmhmm_path, col_types = tmhmm_cols)
lengths <- select(tmhmm_tmp, seqid, len) # For phobius
tmhmm <- process_tmhmm(tmhmm_tmp)


# Phobius
phobius_cols <- cols(
    seqid = col_character(),
    tm = col_integer(),
    sp = col_character(),
    prediction = col_character()
)

phobius <- read_tsv(phobius_path, col_types = phobius_cols) %>%
    process_phobius(lengths)


combined <- bind_rows(
        effectorp,
        apoplastp,
        localizer_eff,
        localizer_plant,
        signalp3_hmm,
        signalp4,
        targetp,
        phobius,
        tmhmm
    ) %>%
    arrange(seqid, analysis, start, end, score)

# Write to stdout
cat(format_tsv(combined))

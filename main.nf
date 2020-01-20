#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include create_db as create_db from './modules/clustering'
include create_db as create_enrich_db from './modules/clustering'
include cluster from './modules/clustering'
include get_cluster_seqs from './modules/clustering'
include split_db as split_msa_db from './modules/misc'
include mafft from './modules/misc'

def help_message() {
    log.info"""
    =================================
    pclust
    =================================

    Usage:

    abaaab

    Mandatory Arguments:
      --seqs              The protein sequences to cluster.
      --db                The proteins to cluster as an MMseqs formatted database.
                          Either `--db` or `--seqs` must be provided to run clustering.
      --enrich_seqs       Sequences to enrich the clusters with for clustering
                          and converting into HHsuite databases.
      --enrich_db         A database to enrich clusters and hhsuite database with.
                          Should be an MMseqs formatted sequence database.
                          If neither `--enrich_seqs` or `--enrich_db` then the seqs
                          will be used for enrichment.

    Options:
      --nomsa             Stop the pipeline after clustering.
      --tree              Predict quick phylogenetic trees for each MSA.
      --noremote_build    Don't build the hhsuite database from the MSAs.
      --noremote          Stop the pipeline after building the hhsuite database,
                          and don't run the HMM-HMM comparisons.
      --clusters          Provide existing clusters as an MMseqs database instead of using
                          the built-in pipeline. Note that the seq database must also be provided.
      --msas              Provide existing MSAs as an MMSeqs MSA database (fasta format).
      --hhself            Provide an existing HHsuite database of clusters to use.
      --hhdata            The path to the HHsuite data folder. If not provided will fetch
                          from the \$HHLIB environment variable.
      --hhpfam            The pfam database formatted as an HHsuite database.
      --hhscop            The scop database formatted as an HHsuite database.
      --hhpdb             The pdb database formatted as an HHsuite database.
      --hhuniref          A uniref (or similar) database formatted as an HHsuite database.
                          Uniclust is a good preformatted option.
      --hhmatches_self
      --hhmatches_pfam
      --hhmatches_scop
      --hhmatches_pdb
      --hhmatches_uniref

    Outputs:

    """.stripIndent()
}


def validate_params(params) {
    def should_run = [:]
    should_run["clustering"] = !params.clusters && !params.msas && !params.hhself
    should_run["msa"] = (should_run["clustering"] || !params.msas) && !params.hhself && !params.nomsa
    should_run["tree"] = (should_run["msa"] || params.msas) && params.tree
    should_run["remote_build"] = (should_run["msa"] || !params.hhself) && !params.noremote_build
    should_run["remote"] = (should_run["remote"] || params.hhself) && !params.noremote

    def error = false

    if ( should_run["clustering"] && !(params.db || params.seqs) ) {
        log.error "The clustering stage " +
                  "requires either '--seqs' or '--db' to be provided."
        error = true
    }

    if ( (should_run["clustering"] || should_run["remote_build"]) && \
         !(params.enrich_db || params.enrich_seqs) ) {
        log.error "An enrichment database is required for the clustering and " +
                  "hhsuite database construction stages. Please specify either " +
                  "'--enrich_seqs' or '--enrich_db'."
        error = true
    }

    if (error) {
        exit 1
    }

    return should_run
}


def param_unexpected_error() {
    log.error "We encountered an error while validating input arguments that " +
        "should be possible. Please raise an issue on github or contact the " +
        "authors."
    exit 1
}


def get_db(db_path, seqs_path, default_db) {

    if ( db_path ) {
        db = Channel.value( file(db_path, checkIfExists: true) )
    } else if ( seqs_path ) {
        seqs = file(seqs_path, checkIfExists: true)
        db = create_db(seqs)
    } else if ( default_db ) {
        db = default_db
    } else {
        param_unexpected_error()
    }
}


def get_file(filepath) {
    if ( filepath ) {
        handle = file(filepath, checkIfExists: true)
    } else {
        param_unexpected_error()
    }

    return handle
}


workflow {

    main:

    if ( params.help ) {
        help_message()
        exit 0
    }

    def should_run = validate_params(params)

    if (should_run["clustering"] || should_run["remote_build"]) {
        if ( params.enrich_db ) {
            enrich_seq_db = Channel.value( file(params.enrich_db, checkIfExists: true) )
        } else if ( params.enrich_seqs ) {
            enrich_seqs = file(params.enrich_seqs, checkIfExists: true)
            enrich_seq_db = create_enrich_db(enrich_seqs)
        } else {
            param_unexpected_error()
        }
    }

    if (should_run["clustering"] || should_run["msa"]) {
        if ( params.db ) {
            seq_db = Channel.value( file(params.db, checkIfExists: true) )
        } else if ( params.seqs ) {
            seqs = file(params.seqs, checkIfExists: true)
            seq_db = create_db(seqs)
        } else {
            param_unexpected_error()
        }
    }

    if (should_run["clustering"]) {
        (cc, cc_tsv, cc_representative, cc_stats,
         pc, pc_tsv, pc_representative, pc_stats) = cluster(seq_db, true, enrich_seq_db)

    } else if ( params.clusters ) {
        pc = file(params.clusters, checkIfExists: true)

        // This junk prevents errors while trying to publish.
        cc = Channel.empty()
        cc_tsv = Channel.empty()
        cc_representative = Channel.empty()
        cc_stats = Channel.empty()
        pc_tsv = Channel.empty()
        pc_representative = Channel.empty()
        pc_stats = Channel.empty()
    }

    if (should_run["msa"]) {
        pc_seqs = get_cluster_seqs(seqs, pc, "profile")
        split_pc_seqs = split_msa_db(20000, pc_seqs).flatten()
        msas = mafft(split_pc_seqs)
    } else {
        pc_seqs = Channel.empty()
        msas = Channel.empty()
    }


    publish:
    cc to: "${params.outdir}/clusters"
    cc_tsv to: "${params.outdir}/clusters"
    cc_representative to: "${params.outdir}/clusters"
    cc_stats to: "${params.outdir}/clusters"
    pc to: "${params.outdir}/clusters"
    pc_tsv to: "${params.outdir}/clusters"
    pc_representative to: "${params.outdir}/clusters"
    pc_stats to: "${params.outdir}/clusters"
    pc_seqs to: "${params.outdir}/clusters"
}

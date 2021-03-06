#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include create_db as create_db from './modules/clustering'
include create_db as create_enrich_db from './modules/clustering'
include cluster from './modules/clustering'
include get_cluster_seqs from './modules/clustering'
include enrich_msas from './modules/clustering'

include check_user_db as check_user_seq_db from './modules/misc'
include check_user_db as check_user_enrich_db from './modules/misc'
include check_user_db as check_user_clusters_db from './modules/misc'
include check_user_db as check_user_msas_db from './modules/misc'
include split_db as split_pc_seqs_db from './modules/misc'
include split_db as split_msa_db from './modules/misc'
include combine_split_dbs as combine_split_msa_dbs from './modules/misc'
include combine_split_dbs as combine_split_tree_dbs from './modules/misc'
include decipher from './modules/misc'

include create_and_sort_hhdb from './modules/hmm'
include copy_hhdata from './modules/hmm'
include split_hhm_databases from './modules/hmm'
include search_seqs_vs_hmms from './modules/hmm'
include search_hmms_vs_hmms as search_self_hmms_vs_hmms from './modules/hmm'
include search_hmms_vs_hmms as search_pfam_hmms_vs_hmms from './modules/hmm'
include search_hmms_vs_hmms as search_scop_hmms_vs_hmms from './modules/hmm'
include search_hmms_vs_hmms as search_pdb_hmms_vs_hmms from './modules/hmm'
include subset_hhdb_by_matches from './modules/hmm'
include combine_hhsuite_results as combine_hmm_matches from './modules/hmm'


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

      --hhmatches_self
      --hhmatches_pfam
      --hhmatches_scop
      --hhmatches_pdb

    Outputs:

    """.stripIndent()
}


def validate_params(params) {
    def should_run = [:]
    should_run["clustering"] = !params.clusters && !params.msas && !params.hhself
    should_run["msa"] = (should_run["clustering"] || !params.msas) && !params.hhself && !params.nomsa
    should_run["tree"] = (should_run["msa"] || params.msas) && params.tree
    should_run["remote_build"] = (should_run["msa"] || !params.hhself) && !params.noremote_build
    should_run["remote"] = (should_run["remote_build"] || params.hhself) && !params.noremote
    should_run["subset_remote"] = (should_run["remote"] && params.hhself_subset_seqs)

    should_run["remote_self"] = (should_run["remote"] && !params.hhmatches_self)
    should_run["remote_pfam"] = (should_run["remote"] && params.hhpfam && !params.hhmatches_pfam)
    should_run["remote_scop"] = (should_run["remote"] && params.hhscop && !params.hhmatches_scop)
    should_run["remote_pdb"] = (should_run["remote"] && params.hhpdb && !params.hhmatches_pdb)

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
        seqs = Channel.value( file(seqs_path, checkIfExists: true) )
        db = create_db(seqs)
    } else if ( default_db ) {
        db = default_db
    } else {
        param_unexpected_error()
    }
}


def get_file(filepath) {
    if ( filepath ) {
        handle = Channel.value( file(filepath, checkIfExists: true) )
    } else {
        // The expectation is that you would check that filepath is not false
        // before you use this.
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

    if (should_run["clustering"] || should_run["msa"]) {
        if ( params.db ) {
            seq_db = check_user_seq_db(get_file(params.db))
        } else if ( params.seqs ) {
            seqs = get_file(params.seqs)
            seq_db = create_db(seqs)
        } else {
            param_unexpected_error()
        }
    }

    if (should_run["clustering"] || should_run["remote_build"]) {
        if ( params.enrich_db ) {
            enrich_seq_db = check_user_enrich_db(get_file(params.enrich_db))
        } else if ( params.enrich_seqs ) {
            enrich_seqs = get_file(params.enrich_seqs)
            enrich_seq_db = create_enrich_db(enrich_seqs)
        } else {
            param_unexpected_error()
        }
    }

    if (params.hhdata) {
        hhdata = Channel.value(get_file(params.hhdata))
    } else if (should_run["remote_build"]) {
        hhdata = copy_hhdata()
    } else {
        hhdata = Channel.empty()
    }


    if (should_run["clustering"]) {
        (cc, cc_tsv, cc_representative, cc_stats,
         pc, pc_tsv, pc_representative, pc_stats) = cluster(seq_db, true, enrich_seq_db)
    } else {
        if ( params.clusters ) {
            pc = check_user_clusters_db(get_file(params.clusters))
        } else {
            pc = Channel.empty()
        }

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
        pc_seqs = get_cluster_seqs(seq_db, pc, "profile")
        split_pc_seqs = split_pc_seqs_db(20000, "protein", pc_seqs)
        split_msas = decipher(split_pc_seqs.flatten())
        msas = combine_split_msa_dbs("msas", split_msas.collect())

    } else {
        if ( params.msas ) {
            msas = check_user_msas_db(get_file(params.msas))
            split_msas = split_msa_db(20000, "protein", msas).chunks.flatten()
        } else {
            msas = Channel.empty()
            split_msas = Channel.empty()
        }

        pc_seqs = Channel.empty()

    }

    if (should_run["tree"]) {
        split_trees = fasttree(split_msas)
        trees = combine_split_tree_dbs("trees", split_trees.collect())
    } else {
        split_trees = Channel.empty()
        trees = Channel.empty()
    }

    if (should_run["remote_build"]) {
        enriched_msas = enrich_msas(split_msas, enrich_seq_db)
        hhself = create_and_sort_hhdb(enriched_msas, hhdata)
    } else {
        hhself = Channel.empty()
        enriched_msas = Channel.empty()
    }

    if (should_run["subset_remote"]) {
        hhself_subset_seqs = get_file(params.hhself_subset_seqs)
        subset_hhself_matches = search_seqs_vs_hmms(
            2,
            0.01,
            20,
            20000,
            10,
            "hhself_subset",
            hhself_subset_seqs,
            hhself
        )
        subset_hhself = subset_hhdb_by_matches(subset_hhself_matches, hhself)

    } else {
        subset_hhself_matches = Channel.empty()
        subset_hhself = hhself
    }

    if (should_run["remote"]) {
        split_hhm = split_hhm_databases(20000, subset_hhself)

        if (should_run["remote_self"]) {
            split_self_hmm_matches = search_self_hmms_vs_hmms(
                1000,
                "split_self",
                split_hhm,
                subset_hhself
            )
            .map { ["hhself", it] }

        } else if (params.hhmatches_self) {
            split_self_hmm_matches = Channel
                .of( ["hhself", get_file(params.hhmatches_self)] )

        } else {
            param_unexpected_error()
        }

        if (should_run["remote_pfam"]) {
            pfam = Channel.value(get_file(params.hhpfam))

            split_pfam_hmm_matches = search_pfam_hmms_vs_hmms(
                1000,
                "split_pfam",
                split_hhm,
                pfam
            )
            .map(["hhpfam", it])

        } else if (params.hhmatches_pfam) {
            split_pfam_hmm_matches = Channel
                .of( ["hhpfam", get_file(params.hhmatches_pfam)] )

        } else {
            split_pfam_hmm_matches = Channel.empty()

        }

        if (should_run["remote_scop"]) {
            scop = Channel.value(get_file(params.hhscop))

            split_scop_hmm_matches = search_scop_hmms_vs_hmms(
                1000,
                "split_scop",
                split_hhm,
                scop
            )
            .map(["hhscop", it])

        } else if (params.hhmatches_scop) {
            split_scop_hmm_matches = Channel
                .of( ["hhscop", get_file(params.hhmatches_scop)] )

        } else {
            split_scop_hmm_matches = Channel.empty()

        }


        if (should_run["remote_pdb"]) {
            pdb = Channel.value(get_file(params.hhpdb))

            split_pdb_hmm_matches = search_pdb_hmms_vs_hmms(
                1000,
                "split_pdb",
                split_hhm,
                pdb
            )
            .map(["hhpdb", it])

        } else if (params.hhmatches_pdb) {
            split_pdb_hmm_matches = Channel
                .of( ["hhpdb", get_file(params.hhmatches_pdb)] )

        } else {
            split_pdb_hmm_matches = Channel.empty()

        }

        (hmm_matches, hmm_matches_tsv) = combine_hmm_matches(
            split_self_hmm_matches.mix(
                split_pfam_hmm_matches,
                split_scop_hmm_matches,
                split_pdb_hmm_matches
            )
        )

    } else {
        split_hhm = Channel.empty()
        self_hmm_matches = Channel.empty()
        hmm_matches = Channel.empty()
        hmm_matches_tsv = Channel.empty()
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
    msas to: "${params.outdir}/msas"
    trees to: "${params.outdir}/msas"
    enriched_msas to: "${params.outdir}/msas"
    hhself to: "${params.outdir}/hhself"
    subset_hhself to: "${params.outdir}/hhself"
    subset_hhself_matches to: "${params.outdir}/hhself"
    hmm_matches to: "${params.outdir}/remote"
    hmm_matches_tsv to: "${params.outdir}/remote"
}

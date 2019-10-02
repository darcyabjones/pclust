#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
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


if (params.help){
    helpMessage()
    exit 0
}


params.seqs = false
params.db = false

params.enrich_db = false
params.enrich_seqs = false

params.enrich_cluster = false

// Options to skip steps, using provided datasets instead
// Skip the clustering.
params.clusters = false
// Skip the multiple sequence alignment.
params.msas = false
// Skip the hhsuite database construction.
params.hhself = false

// Options to stop analyses at points.
params.nomsa = false
params.tree = false
params.noremote_build = false
params.noremote = false

// Remote homology analyses options.
params.hhdata = false

// Provided hhsuite formatted databases to search against.
// These should be in a directory, and the db name prefix should be 'pfam',
// 'scop', 'pdb', or 'uniref'.
// e.g. db/pfam_cs219.ff{data,index} db/pfam_a3m.ff{data,index} etc...
//
// I recognise that having to rename things is painful, but there are too many
// files for each database to provide options for them individually.
params.hhpfam = false
params.hhscop = false
params.hhpdb = false
params.hhuniref = false

// Use these hhblits results instead of running the analyses.
params.hhmatches_self = false
params.hhmatches_pfam = false
params.hhmatches_scop = false
params.hhmatches_pdb = false
params.hhmatches_uniref = false

params.split_size = 20000

def run_clustering = !params.clusters && !params.msas && !params.hhself
def run_msa = (run_clustering || !params.msas) && !params.hhself && !params.nomsa
def run_tree = (run_msa || params.msas) && params.tree
def run_remote_build = (run_msa || !params.hhself) && !params.noremote_build
def run_remote = (run_remote_build || params.hhself) && !params.noremote


//
// STEP 0 - Input validation.
//

if ( params.db ) {

    seqdb = Channel
        .fromPath(
            params.db,
            type: 'dir',
            checkIfExists: true,
            glob: false
        )
        .first()

} else if ( params.seqs ) {

    proteins = Channel
        .fromPath(
            params.seqs,
            type: 'file',
            checkIfExists: true,
            glob: false
        )
        .first()

    /*
     * Create the mmseqs2 sequence database
     */
    process createSequenceDB {

        label 'mmseqs'
        label "small_task"

        publishDir "${params.outdir}/sequences"

        when:
        run_clustering || run_remote_build

        input:
        file "seqs.fasta" from proteins

        output:
        file "seqdb" into seqdb

        script:
        """
        mkdir -p "seqdb"
        mmseqs createdb "seqs.fasta" "seqdb/db" --max-seq-len 14000
        """
    }

} else if ( run_clustering ) {

    log.error "The clustering stage " +
              "requires either '--seqs' or '--db' to be provided."
    exit 1

} else {

    seqdb = Channel.empty()

}


if ( params.enrich_db ) {

    enrichdb = Channel
        .fromPath(
            params.enrich_db,
            type: 'dir',
            checkIfExists: true,
            glob: false
        )
        .first()

} else if ( params.enrich_seqs ) {

    enrichSeqs = Channel
        .fromPath(
            params.enrich_seqs,
            type: 'file',
            checkIfExists: true,
            glob: false
        )
        .first()


    process createEnrichSeqsDB {

        label 'mmseqs'
        label "small_task"

        publishDir "${params.outdir}/sequences"

        when:
        run_clustering || run_remote_build

        input:
        file "seqs.fasta" from enrichSeqs

        output:
        file "enrich_db" into enrichdb

        script:
        """
        mkdir -p "enrich_db"
        mmseqs createdb "seqs.fasta" "enrich_db/db" --max-seq-len 14000
        """
    }

} else if ( params.db || params.seqs ) {

    enrichdb = seqdb

} else {

    log.error "An enrichment database is required for the clustering and " +
              "hhsuite database construction stages. Please specify either " +
              "'--enrich_seqs' or '--enrich_db'."
    exit 1

}


if ( params.hhdata ) {

    hhdata = Channel
        .fromPath(
            params.hhdata,
            type: 'dir',
            checkIfExists: true,
            glob: false
        )
        .first()

} else if ( run_remote_build ) {

    process getHHData {

        label "hhsuite"
        label "small_task"

        when:
        run_remote_build

        output:
        file "hhdata" into hhdata

        script:
        """
        cp -r \${HHLIB}/data hhdata
        """
    }

} else {

    hhdata = Channel.empty()

}


if ( params.hhpfam ) {

    pfamdb = Channel
        .fromPath(
            params.hhpfam,
            type: 'dir',
            checkIfExists: true,
            glob: false
        )
        .first()

} else {

    pfamdb = Channel.empty()

}


if ( params.hhscop ) {

    scopdb = Channel
        .fromPath(
            params.hhscop,
            type: 'dir',
            checkIfExists: true,
            glob: false
        )
        .first()

} else {

    scopdb = Channel.empty()

}


if ( params.hhpdb ) {

    pdbdb = Channel
        .fromPath(
            params.hhpdb,
            type: 'dir',
            checkIfExists: true,
            glob: false
        )
        .first()

} else {

    pdbdb = Channel.empty()

}


if ( params.hhuniref ) {

    unirefdb = Channel
        .fromPath( params.hhuniref, type: 'dir', checkIfExists: true, glob: false )
        .first()

} else {

    unirefdb = Channel.empty()

}


//
// STEP 1. First pass clustering
//


/*
 * Perform the first pass clustering using basic mmseqs workflow.
 */
process clusterCascade {

    label 'mmseqs'
    label "big_task"

    when:
    run_clustering

    input:
    file "seq" from seqdb

    output:
    file "cascade" into cascadeClusters

    script:
    """
    mkdir -p "cascade"
    mkdir -p "tmp"

    mmseqs cluster \
      "seq/db" \
      "cascade/db" \
      "tmp" \
      --threads "${task.cpus}" \
      --min-seq-id 0.3 \
      -c 0.8 \
      --cov-mode 0 \
      --cluster-steps 3 \
      -s 5 \
      --cluster-mode 0 \
      --db-load-mode 0

    rm -rf -- "tmp"
    """
}


/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractCascadeClusterStats {

    label 'mmseqs'
    label "big_task"

    publishDir "${params.outdir}/clusters"

    when:
    run_clustering

    input:
    file "seqs" from seqdb
    file "cascade" from cascadeClusters

    output:
    file "cascade.tsv" into cascadeClustersTSV
    file "cascade_rep.fasta" into cascadeClustersRepFasta
    file "cascade_stats.tsv" into cascadeClustersStats

    script:
    SEQS = "seqs"
    CLUSTERS = "cascade"
    NCPUS = task.cpus
    template "mmseqs_cluster_stats.sh"
}


//
// STEP 2 - second pass clustering using enriched PSSMs.
//

/*
 * Convert the cascade clusters into PSSMs
 */
process createProfile {

    label "mmseqs"
    label "big_task"

    when:
    run_clustering

    input:
    file "seqs" from seqdb
    file "clusters" from cascadeClusters

    output:
    file "profile" into cascadeProfile

    script:
    """
    mkdir -p "profile"

    mmseqs result2profile \
      "seqs/db" \
      "seqs/db" \
      "clusters/db" \
      "profile/db" \
      --threads "${task.cpus}"
    """
}


/*
 * Enrich the profile by searching a database.
 */
process enrichProfile {

    label "mmseqs"
    label "big_task"

    when:
    run_clustering && ( params.enrich_db || params.enrich_seqs ) && params.enrich_cluster

    input:
    file "profile" from cascadeProfile
    file "enrich_seqs" from enrichdb
    file "seqs" from seqdb

    output:
    file "enrich_matches" into enrichedSearchResults

    script:
    """
    mkdir -p "tmp"
    mkdir -p "enrich_matches"

    mmseqs search \
      "profile/db" \
      "enrich_seqs/db" \
      "enrich_matches/db" \
      "tmp" \
      --threads "${task.cpus}" \
      --max-seqs 300 \
      -e 0.00001 \
      -s 5 \
      --rescore-mode 1 \
      --db-load-mode 0 \
      --split 0

    rm -rf -- "tmp"
    """
}


/*
 * Convert search results into an enriched profile.
 */
process createEnrichedProfile {

    label "mmseqs"
    label "big_task"

    when:
    run_clustering && ( params.enrich_db || params.enrich_seqs ) && params.enrich_cluster

    input:
    file "input_profile" from cascadeProfile
    file "enrich_seqs" from enrichdb
    file "enrich_matches" from enrichedSearchResults

    output:
    file "enriched_profile" into enrichedProfile

    script:
    """
    mkdir -p "enriched_profile"

    mmseqs result2profile \
      "input_profile/db" \
      "enrich_seqs/db" \
      "enrich_matches/db" \
      "enriched_profile/db" \
      --threads "${task.cpus}"
    """
}


if ( (params.enrich_db || params.enrich_seqs) && params.enrich_cluster ) {

    enrichedProfile.into { profile4CluSearch; profile4Clu }

} else {

    cascadeProfile.into { profile4CluSearch; profile4Clu }

}


/*
 * Search the cluster profiles against the cluster consensus sequences.
 */
process clusterProfileSearch {

    label 'mmseqs'
    label "big_task"

    when:
    run_clustering

    input:
    file "input_profile" from profile4CluSearch

    output:
    file "profile_matches" into profileClusterSearchResults

    script:
    """
    mkdir "tmp"
    mkdir "profile_matches"

    mmseqs search \
      "input_profile/db" \
      "input_profile/db_consensus" \
      "profile_matches/db" \
      "tmp" \
      --threads "${task.cpus}" \
      --max-seqs 100 \
      -c 0.8 \
      --min-seq-id 0.3 \
      --cov-mode 1 \
      -s 6.0 \
      -e 0.00001 \
      --db-load-mode 0 \
      --split 0 \
      --add-self-matches

    rm -rf -- "tmp"
    """
}


/*
 * Cluster the cluster profiles based on the profile/consensus search.
 */
process clusterProfile {

    label 'mmseqs'
    label "big_task"

    when:
    run_clustering

    input:
    file "input_profile" from profile4Clu
    file "profile_matches" from profileClusterSearchResults

    output:
    file "profile" into profileClusters

    script:
    """
    mkdir "profile"
    mmseqs clust \
      "input_profile/db" \
      "profile_matches/db" \
      "profile/db" \
      --threads "${task.cpus}" \
      --cluster-mode 2
    """

}


/*
 * Merge the cascade and profile clustering results to get final clusters.
 */
process mergeClusters {

    label "mmseqs"
    label "small_task"

    publishDir "${params.outdir}/clusters"

    when:
    run_clustering

    input:
    file "seq" from seqdb
    file "cascade" from cascadeClusters
    file "profile_tmp" from profileClusters

    output:
    file "profile" into mergedClusters
    file "profile_seqs" into mergedClustersSeqs

    """
    mkdir "profile"

    mmseqs mergeclusters \
      "seq/db" \
      "profile/db" \
      "cascade/db" \
      "profile_tmp/db"

    mkdir "profile_seqs"
    mmseqs createseqfiledb \
      "seq/db" \
      "profile/db" \
      "profile_seqs/db"

    ORIG="\${PWD}"
    cd profile_seqs
    ln -s db db.ffdata
    ln -s db.dbtype db.ffdata.dbtype
    ln -s db.index db.ffindex
    cd \${ORIG}
    """
}


/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractProfileClusterStats {

    label 'mmseqs'
    label "big_task"

    publishDir "${params.outdir}/clusters"

    when:
    run_clustering

    input:
    file "profile" from mergedClusters
    file "seq" from seqdb

    output:
    file "profile.tsv" into profileClustersTSV
    file "profile_rep.fasta" into profileClustersRepFasta
    file "profile_stats.tsv" into profileClustersStats

    script:
    SEQS = "seq"
    CLUSTERS = "profile"
    NCPUS = task.cpus
    template "mmseqs_cluster_stats.sh"
}


if ( params.clusters ) {

    clusters = Channel
        .fromPath(
            params.clusters,
            type: 'dir',
            checkIfExists: true,
            glob: false
        )
        .first()

} else {

    clusters = mergedClustersSeqs

}


//
// STEP 3. Get multiple sequence alignments of each cluster.
//

/*
 * Get fasta sequence databases from the cluster database and
 * and split the database into many parts to allow checkpointing.
 */
process splitClusterSeqs {

    label "mmseqs"
    label "big_task"

    when:
    run_msa

    input:
    file "clusters" from clusters

    output:
    file "split_seqs_*" into splitClusters mode flatten

    script:
    """
    TARGET_CLUSTER_SIZE="${params.split_size}"
    NCLUSTERS=\$(wc -l < "clusters/db.index")
    NSPLITS=\$(( (\${NCLUSTERS} + \${TARGET_CLUSTER_SIZE} + 1) / \${TARGET_CLUSTER_SIZE} ))

    mmseqs splitdb "clusters/db" "tmp_split_seqs" --split \${NSPLITS}

    for f in tmp_split_seqs_*.index
    do
        BASENAME="\${f%.index}"
        DIRNAME="\${BASENAME#tmp_}"

        mkdir "\${DIRNAME}"
        mv "\${f}" "\${DIRNAME}/db.index"
        mv "\${BASENAME}" "\${DIRNAME}/db"
        mv "\${BASENAME}.dbtype" "\${DIRNAME}/db.dbtype"
    done
    """
}


/*
 * Construct MSA from the fasta databases.
 */
process mafftMSA {

    label "mafft"
    label "big_task"

    when:
    run_msa

    input:
    file "inmsas" from splitClusters

    output:
    file "msas" into splitMSAs

    script:
    """
    mkdir -p "msas"
    mpirun -np "${task.cpus}" mmseqs apply \
      "inmsas/db" \
      "msas/db" \
      --threads 1 \
      -- \
      run_mafft.sh
    """
}

splitMSAs.into {
    splitMSAs4CombineSplitMSAs;
    splitMSAs4DecideIfUser;
}


/*
 * Collect the split MSA databases into a single MSA database.
 */
process combineSplitMSAs {

    label "mmseqs"
    label "small_task"

    publishDir "${params.outdir}"

    when:
    run_msa

    input:
    file "split_msas_*" from splitMSAs4CombineSplitMSAs.collect()

    output:
    file "msas" into combinedMSAs

    script:
    """
    for db in split_msas_*
    do
      if [ -e "msas" ]
      then
        mkdir "reduce"
        mmseqs concatdbs "msas/db" "\${db}/db" "reduce/db" --preserve-keys
        rm -rf -- "msas"
        mv "reduce" "msas"
      else
        cp -r -L "\${db}" "msas"
      fi
    done

    ORIG=\${PWD}
    cd msas
    ln -s "db" "db.ffdata"
    ln -s "db.dbtype" "db.ffdata.dbtype"
    ln -s "db.index" "db.ffindex"
    cd "\${ORIG}"
    """
}


if ( params.msas ) {

    userMsa = Channel
        .fromPath(
            params.msas,
            type: 'dir',
            checkIfExists: true,
            glob: false
        )
        .first()


    process processUserMSA {

        label "mmseqs"
        label "small_task"

        input:
        file "inmsas" from userMsa

        output:
        file "msas" into msas
        file "split_msa_*" into splitUserMSAs mode flatten

        script:
        """
        cp -r -L inmsas msas

        ORIG=\${PWD}

        if [ ! -e msas/db ] && [ -e msas/db.ffdata ]
        then
            cd msas
            ln -s db.ffdata db
            cd "\${ORIG}"
        elif [ ! -e msas/db.ffdata ] && [ -e msas/db ]
        then
            cd msas
            ln -s db db.ffdata
            cd "\${ORIG}"
        elif [ -e msas/db ] && [ -e msas/db.ffdata ]
        then
            true
        else
            echo "An msa was specified but there was no db or db.ffdata file!" 1>&2
            exit 1
        fi


        if [ ! -e msas/db.index ] && [ -e msas/db.ffindex ]
        then
            cd msas
            ln -s db.ffindex db.index
            cd "\${ORIG}"
        elif [ ! -e msas/db.ffindex ] && [ -e msas/db.index ]
        then
            cd msas
            ln -s db.index db.ffindex
            cd "\${ORIG}"
        elif [ -e msas/db ] && [ -e msas/db.ffdata ]
        then
            true
        else
            echo "An msa was specified but there was no db.index or db.ffindex file!" 1>&2
            exit 1
        fi


        # split the db
        TARGET_CLUSTER_SIZE="${params.split_size}"
        NCLUSTERS=\$(wc -l < "msas/db.index")
        NSPLITS=\$(( (\${NCLUSTERS} + \${TARGET_CLUSTER_SIZE} + 1) / \${TARGET_CLUSTER_SIZE} ))

        mmseqs splitdb "msas/db" "tmp_split_msa" --split \${NSPLITS}

        for f in tmp_split_msa_*.index
        do
            BASENAME="\${f%.index}"
            DIRNAME="\${BASENAME#tmp_}"

            mkdir "\${DIRNAME}"
            mv "\${f}" "\${DIRNAME}/db.index"
            mv "\${BASENAME}" "\${DIRNAME}/db"

            cd "\${DIRNAME}"
            ln -s "db.index" "db.ffindex"
            ln -s "db" "db.ffdata"
            cd "\${ORIG}"
        done
        """
    }

    splitUserMSAs.into {
        splitMSAs4EnrichMSA;
        splitMSAs4ConstructTree;
    }

} else {

    msas = combinedMSAs
    splitMSAs4DecideIfUser.into {
        splitMSAs4EnrichMSA;
        splitMSAs4ConstructTree;
    }

}


//
// STEP 4 - Get a tree for each cluster MSA
//


/*
 * Construct an ML tree for each cluster.
 */
process constructTree {

    label "fasttree"
    label "big_task"

    when:
    run_tree

    input:
    file "msas" from splitMSAs4ConstructTree

    output:
    file "trees" into splitTrees

    script:
    """
    OMP_NUM_THREADS=1

    mkdir -p "trees"
    mpirun -np "${task.cpus}" mmseqs apply \
      "msas/db" \
      "trees/db" \
      --threads 1 \
      -- \
      FastTree -quiet -nopr
    """
}


/*
 * Collect the split tree databases into a single tree database.
 */
process combineSplitTrees {

    label "mmseqs"
    label "small_task"

    publishDir "${params.outdir}"

    when:
    run_tree

    input:
    file "split_trees_*" from splitTrees.collect()

    output:
    file "trees" into combinedTrees

    script:
    """
    for db in split_trees_*
    do
      if [ -e "trees" ]
      then
        mkdir "reduce"
        mmseqs concatdbs "trees/db" "\${db}/db" "reduce/db" --preserve-keys
        rm -rf -- "trees"
        mv "reduce" "trees"
      else
        cp -r -L "\${db}" "trees"
      fi
    done

    ORIG="\${PWD}"

    cd trees
    ln -s "db" "db.ffdata"
    ln -s "db.dbtype" "db.ffdata.dbtype"
    ln -s "db.index" "db.ffindex"

    cd "\${ORIG}"
    """
}


//
// STEP 5 - Convert the MSAs into an HHBlits database.
//

/*
 * Enrich the MSAs by searching against a database.
 * This will give extra information that might not be there
 * because we have quite strict coverage requirements.
 * It will help matching domains later on.
 *
 * We use MMSeqs because it is MUCH faster than the HHBlits enrichment
 * method.
 */
process enrichMSA {

    label "mmseqs"
    label "big_task"

    when:
    run_remote_build

    input:
    file "msas" from splitMSAs4EnrichMSA
    file "enrich" from enrichdb

    output:
    file "enriched_msas" into enrichedMSAs

    script:
    """
    mkdir "profile"
    mmseqs msa2profile \
      "msas/db" \
      "profile/db" \
      --match-mode 1 \
      --match-ratio 1

    mkdir "search" "tmp"
    mmseqs search \
      "profile/db" \
      "enrich/db" \
      "search/db" \
      "tmp" \
      --threads "${task.cpus}" \
      -s 6 \
      -a \
      -e 0.00001 \
      --db-load-mode 0 \
      --split 0 \
      --split-mode 1


    mkdir "search_msas"
    mmseqs result2msa \
      "profile/db" \
      "enrich/db" \
      "search/db" \
      "search_msas/db" \
      --skip-query 1

    mkdir "enriched_msas"
    ffdb join_concat \
      -d "enriched_msas/db" \
      -i "enriched_msas/db.index" \
      "msas/db" "search_msas/db" \
      "msas/db.index" "search_msas/db.index"

    # Would be good to sort by num msa columns to help balance load.
    # Not necessary short-term, order seems like that anyway.

    cp -L "search_msas/db.dbtype" "enriched_msas/db.dbtype"

    ORIG="\${PWD}"
    cd "enriched_msas"
    ln -s "db" "db.ffdata"
    ln -s "db.dbtype" "db.ffdata.dbtype"
    ln -s "db.index" "db.ffindex"

    cd "\${ORIG}"

    rm -rf -- "profile" "search" "search_msas" "tmp"
    """
}


/*
 * Create HHsuite databases for each MSA chunk.
 */
process fasToHHDB {

    label "hhsuite"
    label "big_task"

    when:
    run_remote_build

    input:
    file "enriched_msas" from enrichedMSAs
    file "hhdata" from hhdata

    output:
    file "hhdb" into splitHHDB

    script:
    """
    mkdir -p hhdb

    # This should probably get moved to the enrich MSA task for next run.
    mpirun -np "${task.cpus}" ffindex_apply_mpi \
        enriched_msas/db.ff{data,index} \
        -i "hhdb/db_fasta.ffindex" \
        -d "hhdb/db_fasta.ffdata" \
        -- \
        tidy_joined_alignments.py

    cp -L "enriched_msas/db.ffdata.dbtype" "hhdb/db_fasta.ffdata.dbtype"

    mpirun -np "${task.cpus}" ffindex_apply_mpi \
        hhdb/db_fasta.ff{data,index} \
        -i "hhdb/db_a3m.ffindex" \
        -d "hhdb/db_a3m.ffdata" \
        -- \
        run_fas_to_a3m.sh

    mpirun -np "${task.cpus}" cstranslate_mpi \
        -i "hhdb/db_a3m" \
        -o "hhdb/db_cs219" \
        -x 0.3 \
        -c 4 \
        -b \
        -I a3m \
        -A "hhdata/cs219.lib" \
        -D "hhdata/context_data.lib"

    mpirun -np "${task.cpus}" ffindex_apply_mpi \
        hhdb/db_a3m.ff{data,index} \
        -i "hhdb/db_hhm.ffindex" \
        -d "hhdb/db_hhm.ffdata" \
        -- \
        hhmake -i stdin -o stdout -v 0
    """
}


splitHHDB.into {
    splitHHDB4CombineSplitHHDBs;
    splitHHDB4DecideIfUser;
}


/*
 * Combine the HHsuite database chunks into a final database.
 * We can use this database to search against.
 */
process combineSplitHHDBs {

    label "hhsuite"
    label "small_task"

    when:
    run_remote_build

    input:
    file "split_hhdata_*" from splitHHDB4CombineSplitHHDBs.collect()

    output:
    set val("cs219"),
        file("order.txt"),
        file("db_cs219.ffdata"),
        file("db_cs219.ffindex") into hhselfCS219

    set val("fasta"),
        file("order.txt"),
        file("db_fasta.ffdata"),
        file("db_fasta.ffindex") into hhselfFasta

    set val("a3m"),
        file("order.txt"),
        file("db_a3m.ffdata"),
        file("db_a3m.ffindex") into hhselfA3m

    set val("hhm"),
        file("order.txt"),
        file("db_hhm.ffdata"),
        file("db_hhm.ffindex") into hhselfHHM

    script:
    """
    cat split_hhdata_*/db_cs219.ffindex \
    | sort -k3,3n \
    | cut -f1 \
    > order.txt

    ffdb combine \
      -d "db_cs219.ffdata" \
      -i "db_cs219.ffindex" \
      split_hhdata_*/db_cs219.ff{data,index}

    ffdb combine \
      -d "db_fasta.ffdata" \
      -i "db_fasta.ffindex" \
      split_hhdata_*/db_fasta.ff{data,index}

    ffdb combine \
      -d "db_a3m.ffdata" \
      -i "db_a3m.ffindex" \
      split_hhdata_*/db_a3m.ff{data,index}

    ffdb combine \
      -d "db_hhm.ffdata" \
      -i "db_hhm.ffindex" \
      split_hhdata_*/db_hhm.ff{data,index}
    """
}


process sortCombinedHHDBs {

    label "hhsuite"
    label "small_task"
    
    tag "${kind}"

    when:
    run_remote_build

    input:
    set val(kind),
        file("order.txt"),
        file("db.ffdata"),
        file("db.ffindex") from hhselfCS219
            .mix( hhselfFasta, hhselfA3m, hhselfHHM )

    output:
    set file("db_${kind}.ffdata"),
        file("db_${kind}.ffindex") into hhselfSorted

    script:
    """
    ffindex_order \
      "order.txt" \
      db.ff{data,index} \
      db_${kind}.ff{data,index}
    """
}


process combineSortedHHDBs {

    label "posix"
    label "small_task"

    publishDir "${params.outdir}"

    when:
    run_remote_build

    input:
    file "*" from hhselfSorted.collect() 

    output:
    file "hhself" into HHDBTMP

    script:
    """
    mkdir hhself
    ln -sf "\${PWD}/"*.ff{data,index} "\${PWD}/hhself"
    """
}


if ( params.hhself ) {

    userHHSelf = Channel
        .fromPath(
            params.hhself,
            type: 'dir',
            checkIfExists: true,
            glob: false
        )
        .first()


    process splitUserHHSelf {

        label "ffdb"
        label "small_task"

        input:
        file "hhself" from userHHSelf

        output:
        file "hhself" into HHDB
        file "split_hhself_*" into splitUserHHDB

        script:
        """
        # split the db
        TARGET_CLUSTER_SIZE="${params.split_size}"

        ffdb split \
          --size "\${TARGET_CLUSTER_SIZE}" \
          --basename "tmp_split_hhself_{index}.{ext}" \
          "hhself/db_hhm.ffdata" \
          "hhself/db_hhm.ffindex"

        for f in tmp_split_hhself_*.ffindex
        do
            BASENAME="\${f%.ffindex}"
            DIRNAME="\${BASENAME#tmp_}"

            mkdir "\${DIRNAME}"
            mv "\${f}" "\${DIRNAME}/db_hhm.ffindex"
            mv "\${BASENAME}.ffdata" "\${DIRNAME}/db_hhm.ffdata"
        done
        """

    }

    splitUserHHDB.into {
        splitHHDB4SearchSelf;
        splitHHDB4SearchPfam;
        splitHHDB4SearchScop;
        splitHHDB4SearchPdb;
        splitHHDB4SearchUniref;
    }

} else {

    HHDB = HHDBTMP

    splitHHDB4DecideIfUser.into {
        splitHHDB4SearchSelf;
        splitHHDB4SearchPfam;
        splitHHDB4SearchScop;
        splitHHDB4SearchPdb;
        splitHHDB4SearchUniref;
    }

}


//
// STEP 6 - Compare HMMs to find relationships between clusters.
//


/*
 * Perform all-vs-all matches against self.
 * This is to define relationships between clusters.
 */
process searchSelf {

    label "hhsuite"
    label "big_task"

    when:
    run_remote && !params.hhmatches_self

    input:
    set file("subdb"), file("hhdb") from splitHHDB4SearchSelf
        .combine(HHDB)

    output:
    file "results" into splitSelfMatchesTmp

    script:
    """
    mkdir "results"
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i "subdb/db_hhm" \
      -d "hhdb/db" \
      -o "results/db_hhr" \
      -n 1 \
      -cpu 1 \
      -v 0 \
      -e 0.001 \
      -E 0.001 \
      -z 0 \
      -Z 20000 \
      -b 0 \
      -B 20000 \
      -pre_evalue_thresh 100 \
      -min_prefilter_hits 100 \
      -realign_max 20000
    """
}


if ( params.hhmatches_self ) {

    splitSelfMatches = Channel.fromPath(
        params.hhmatches_self,
        type: 'dir',
        checkIfExists: true,
        glob: true
    )

} else {

    splitSelfMatches = splitSelfMatchesTmp

}


/*
 * Search for hmm matches to the PFAM database.
 */
process searchPfam {

    label "hhsuite"
    label "big_task"

    when:
    run_remote && !params.hhmatches_pfam

    input:
    set file("subdb"), file("pfam") from splitHHDB4SearchPfam
        .combine(pfamdb)

    output:
    file "results" into splitPfamMatchesTmp

    script:
    """
    mkdir "results"
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i "subdb/db_hhm" \
      -d "pfam/pfam" \
      -o "results/db_hhr" \
      -n 1 \
      -cpu 1 \
      -v 0 \
      -e 0.001 \
      -E 0.001 \
      -z 0 \
      -Z 500 \
      -b 0 \
      -B 500 \
      -pre_evalue_thresh 10 \
      -min_prefilter_hits 10 \
      -realign_max 500
    """
}

if ( params.hhmatches_pfam ) {

    splitPfamMatches = Channel.fromPath(
        params.hhmatches_pfam,
        type: 'dir',
        checkIfExists: true,
        glob: true
    )

} else {

    splitPfamMatches = splitPfamMatchesTmp

}


/*
 * Search for hmm matches to the SCOP database.
 */
process searchScop {

    label "hhsuite"
    label "big_task"

    when:
    run_remote && !params.hhmatches_scop

    input:
    set file("subdb"), file("scop") from splitHHDB4SearchScop
        .combine(scopdb)

    output:
    file "results" into splitScopMatchesTmp

    script:
    """
    mkdir "results"
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i "subdb/db_hhm" \
      -d "scop/scop" \
      -o "results/db_hhr" \
      -n 1 \
      -cpu 1 \
      -v 0 \
      -e 0.001 \
      -E 0.001 \
      -z 0 \
      -Z 500 \
      -b 0 \
      -B 500 \
      -pre_evalue_thresh 10 \
      -min_prefilter_hits 10 \
      -realign_max 500
    """
}


if ( params.hhmatches_scop ) {

    splitScopMatches = Channel.fromPath(
        params.hhmatches_scop,
        type: 'dir',
        checkIfExists: true,
        glob: true
    )

} else {

    splitScopMatches = splitScopMatchesTmp

}


/*
 * Search for hmm matches to the PDB database.
 */
process searchPdb {

    label "hhsuite"
    label "big_task"

    when:
    run_remote && !params.hhmatches_pdb

    input:
    set file("subdb"), file("pdb") from splitHHDB4SearchPdb
        .combine(pdbdb)

    output:
    file "results" into splitPdbMatchesTmp

    script:
    """
    mkdir "results"
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i "subdb/db_hhm" \
      -d "pdb/pdb" \
      -o "results/db_hhr" \
      -n 1 \
      -cpu 1 \
      -v 0 \
      -e 0.001 \
      -E 0.001 \
      -z 10 \
      -Z 500 \
      -b 10 \
      -B 500 \
      -pre_evalue_thresh 10 \
      -min_prefilter_hits 100 \
      -realign_max 500
    """
}

if ( params.hhmatches_pdb ) {

    splitPdbMatches = Channel.fromPath(
        params.hhmatches_pdb,
        type: 'dir',
        checkIfExists: true,
        glob: true
    )

} else {

    splitPdbMatches = splitPdbMatchesTmp

}


/*
 * Search for hmm matches to a uniref database e.g. uniclust.
 */
process searchUniref {

    label "hhsuite"
    label "big_task"

    when:
    run_remote && !params.hhmatches_uniref

    input:
    set file("subdb"), file("uniref") from splitHHDB4SearchUniref
        .combine(unirefdb)

    output:
    file "results" into splitUnirefMatchesTmp

    script:
    """
    mkdir "results"
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i "subdb/db_hhm" \
      -d "uniref/uniref" \
      -o "results/db_hhr" \
      -n 1 \
      -cpu 1 \
      -v 0 \
      -e 0.001 \
      -E 0.001 \
      -z 0 \
      -Z 100 \
      -b 0 \
      -B 100 \
      -pre_evalue_thresh 1 \
      -min_prefilter_hits 10 \
      -realign_max 100
    """
}


if ( params.hhmatches_uniref ) {

    splitUnirefMatches = Channel.fromPath(
        params.hhmatches_uniref,
        type: 'dir',
        checkIfExists: true,
        glob: true
    )

} else {

    splitUnirefMatches = splitUnirefMatchesTmp

}


splitSelfMatches
    .collect().map { ds -> ["self", ds] }
    .mix(
        splitPfamMatches.collect().map { ds -> ["pfam", ds] },
        splitScopMatches.collect().map { ds -> ["scop", ds] },
        splitPdbMatches.collect().map { ds -> ["pdb", ds] },
        splitUnirefMatches.collect().map { ds -> ["uniref", ds] }
    )
    .set { splitMatches }


/*
 * Combines the split search results into a single database.
 */
process combineSearchMatches {

    label "ffdb"
    label "small_task"

    publishDir "${params.outdir}/hhmatches"

    tag "${database}"

    input:
    set val(database), file("split_results_*") from splitMatches

    output:
    set val(database), file("${database}") into combinedMatches

    script:
    """
    mkdir -p "${database}"

    ffdb combine \
      -d "${database}/db_hhr.ffdata" \
      -i "${database}/db_hhr.ffindex" \
      split_results_*/db_hhr.ff{data,index}
    """
}

#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    =================================
    pclust/pclust
    =================================

    Usage:

    abaaab

    Mandatory Arguments:
      --seqs              description
      --db
      --enrich_seqs
      --enrich_db

    Options:
      --nomsa
      --nomsa_refine
      --enrich_seqs
      --enrich_db
      --enrich_msa

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

// Options to skip steps, using provided datasets instead
// Skip the clustering.
params.clusters = false
// Skip the multiple sequence alignment.
params.msas = false
// Skip the hhsuite database construction.
params.hhself = false

// Options to stop analyses at points.
params.nomsa = false
params.notree = false
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


def run_clustering = !params.clusters || !params.msas || !params.hhself
def run_msa = (!params.msas || !params.hhself) && !params.nomsa
def run_tree = (params.msas || run_msa) && !params.notree
def run_remote = !params.hhself && !params.noremote

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

        input:
        file "seqs.fasta" from proteins

        output:
        file "seqdb" into seqdb

        script:
        """
        mkdir -p "seqdb"
        mmseqs createdb "seqs.fasta" "seqdb/db" --max-seq-len 14000

        # mkdir -p tmp
        # mmseqs createindex "seqdb/db" tmp --threads "${task.cpus}"
        # rm -rf -- tmp
        """
    }

} else if ( run_clustering && run_msa ) {

    log.error "The clustering and multiple sequence alignment stages " +
              "require either '--seqs' or '--db' to be provided."
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

        input:
        file "seqs.fasta" from enrichSeqs

        output:
        file "enrich_db" into enrichdb

        script:
        """
        mkdir -p "enrich_db"
        mmseqs createdb "seqs.fasta" "enrich_db/db" --max-seq-len 14000

        # mkdir -p tmp
        # mmseqs createindex "enrich_db/db" tmp --threads "${task.cpus}"
        # rm -rf -- tmp
        """
    }

} else if ( param.db || params.seqs ) {

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

} else if ( run_remote ) {

    process getHHData {

        label "hhsuite"
        label "small_task"

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
    mkdir -p tmp

    # mmseqs touchdb "${INDB}/db" --threads "${task.cpus}"

    mmseqs cluster \
      "seq/db" \
      "cascade/db" \
      tmp \
      --threads "${task.cpus}" \
      --min-seq-id 0.0 \
      -c 0.8 \
      --cov-mode 0 \
      --cluster-steps 3 \
      -s 6.5 \
      --cluster-mode 0 \
      --db-load-mode 0

    rm -rf -- tmp
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

    publishDir "${params.outdir}/clusters"

    when:
    run_clustering

    input:
    file "clusters" from cascadeClusters
    file "seqs" from seqdb

    output:
    file "profile" into cascadeProfile

    script:
    QUERY = "seqs"
    TARGET = "seqs"
    RESULTS = "clusters"
    OUTDB = "profile"
    NCPUS = task.cpus
    template "mmseqs_result_to_profile.sh"
}


/*
 * Enrich the profile by searching a database.
 */
process enrichProfile {
    label "mmseqs"
    label "big_task"

    when:
    run_clustering && ( params.enrich_db || params.enrich_seqs )

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
      tmp \
      --threads "${task.cpus}" \
      --max-seqs 300 \
      -e 0.00001 \
      --e-profile 0.01 \
      --start-sens 4.5 \
      --sens-steps 2 \
      -s 7.5 \
      --rescore-mode 1 \
      --db-load-mode 0 \
      --split 0

    rm -rf -- tmp
    """
}


/*
 * Convert search results into an enriched profile.
 */
process createEnrichedProfile {
    label "mmseqs"
    label "big_task"

    when:
    run_clustering && ( params.enrich_db || params.enrich_seqs )

    input:
    file "input_profile" from cascadeProfile
    file "enrich_seqs" from enrichdb
    file "enrich_matches" from enrichedSearchResults

    output:
    file "enriched_profile" into enrichedProfile

    script:
    QUERY = "input_profile"
    TARGET = "enrich_seqs"
    RESULTS = "enrich_matches"
    OUTDB = "enriched_profile"
    NCPUS = task.cpus
    template "mmseqs_result_to_profile.sh"
}


if ( params.enrich_db || params.enrich_seqs ) {

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

    publishDir "${params.outdir}/clusters"

    when:
    run_clustering

    input:
    file "input_profile" from profile4CluSearch

    output:
    file "profile_matches" into profileClusterSearchResults
    file "profile_matches.tsv" into profileClusterSearchResultsTSV

    script:
    """
    mkdir -p tmp
    mkdir profile_matches

    mmseqs search \
      "input_profile/db" \
      "input_profile/db_consensus" \
      "profile_matches/db" \
      "tmp" \
      --threads "${task.cpus}" \
      --max-seqs 100 \
      -c 0.8 \
      --cov-mode 0 \
      --start-sens 4 \
      --sens-steps 3 \
      -s 7.5 \
      -e 0.00001 \
      --e-profile 0.01 \
      --db-load-mode 0 \
      --split 0 \
      --add-self-matches

    mmseqs convertalis \
      "input_profile/db" \
      "input_profile/db_consensus" \
      "profile_matches/db" \
      "profile_matches.tsv" \
      --threads "${task.cpus}" \
      --format-mode 0 \
      --format-output "query,target,evalue,qcov,tcov,gapopen,pident,nident,mismatch,raw,bits,qstart,qend,tstart,tend,qlen,tlen,alnlen"

    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' "profile_matches.tsv"

    rm -rf -- tmp
    """
}


/*
 * Cluster the cluster profiles based on the profile/consensus search.
 */
process clusterProfile {

    label 'mmseqs'
    label "big_task"

    publishDir "${params.outdir}/clusters"

    when:
    run_clustering

    input:
    file "input_profile" from profile4Clu
    file "profile_matches" from profileClusterSearchResults

    output:
    file "profile" into profileClusters

    script:
    """
    mkdir -p "profile"
    mmseqs clust \
      "input_profile/db" \
      "profile_matches/db" \
      "profile/db" \
      --threads "${task.cpus}" \
      --cluster-mode 0
    """

}


/*
 * Merge the cascade and profile clustering results to get final clusters.
 */
process mergeClusters {
    label "mmseqs"
    label "big_task"

    publishDir "${params.outdir}/clusters"

    when:
    run_clustering

    input:
    file "seq" from seqdb
    file "cascade" from cascadeClusters
    file "profile" from profileClusters

    output:
    file "merged" into mergedClusters

    """
    mkdir -p merged

    mmseqs mergeclusters \
      seq/db \
      merged/db \
      cascade/db \
      profile/db
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

    clusters = mergedClusters

}


//
// STEP 3. Get multiple sequence alignments of each cluster.
//

/*
 * Get fasta sequence databases from the cluster database and
 * and split the database into many parts to allow checkpointing.
 */
process clusterSeqdb {
    label "mmseqs"
    label "big_task"

    when:
    run_msa

    input:
    file "clusters" from clusters
    file "seq" from seqdb

    output:
    file "cluster_seqs_*" into splitClusters mode flatten

    script:
    """
    mkdir -p cluster_seqs
    mmseqs createseqfiledb "seq/db" "clusters/db" "cluster_seqs/db"

    TARGET_CLUSTER_SIZE=10000
    NCLUSTERS=\$(wc -l < "cluster_seqs/db.index")
    NSPLITS=\$(( (\${NCLUSTERS} + \${TARGET_CLUSTER_SIZE} + 1) / \${TARGET_CLUSTER_SIZE} ))

    mmseqs splitdb "cluster_seqs/db" "cluster_seqs" --split \${NSPLITS}
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
    set file(db), file(db_type), file(db_index) from splitClusters
        .map { [it.getSimpleName(), it] }
        .groupTuple(by: 0, size: 3)
        .map { bn, files -> files }

    output:
    file "msas_${db.getName()}" into splitMSAs

    script:
    """
    mkdir -p "msas_${db.getName()}"
    mpirun -np "${task.cpus}" mmseqs apply \
      "${db}" \
      "msas_${db.getName()}/db" \
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

    publishDir "${params.outdir}/msas"

    when:
    run_msa

    input:
    file "*" from splitMSAs4CombineSplitMSAs.collect()

    output:
    file "msa" into combinedMSAs

    script:
    """
    for db in \$(find . -name "msas*")
    do
      if [ -e "msa" ]
      then
        mkdir "reduce"
        mmseqs concatdbs "msa/db" "\${db}/db" "reduce/db" --preserve-keys
        rm -rf -- "msa"
        mv "reduce" "msa"
      else
        cp -r -L "\${db}" "msa"
      fi
    done

    ln -s msa/db msa/db.ffdata
    ln -s msa/db.dbtype msa/db.ffdata.dbtype
    ln -s msa/db.index msa/db.ffindex
    """
}


if ( params.msa ) {

    userMsa = Channel
        .fromPath(
            params.msa,
            type: 'dir',
            checkIfExists: true,
            glob: false
        )
        .first()


    process processUserMSA {

        label "mmseqs"
        label "small_task"

        input:
        file "inmsa" from userMsa

        output:
        file "msa" into msa
        file "split_msa_*" splitMSAs4EnrichMSA

        script:
        """
        cp -r -L inmsa msa

        if [ ! -e msa/db && -e msa/db.ffdata ]
        then
            ln -s msa/db.ffdata msa/db
        elif [ ! -e msa/db.ffdata && -e msa/db ]
        then
            ln -s msa/db msa/db.ffdata
        elif [ -e msa/db && -e msa/db.ffdata ]
        then
            true
        else
            echo "An msa was specified but there was no db or db.ffdata file!" 1>&2
            exit 1
        fi


        if [ ! -e msa/db.index && -e msa/db.ffindex ]
        then
            ln -s msa/db.ffindex msa/db.index
        elif [ ! -e msa/db.ffindex && -e msa/db.index ]
        then
            ln -s msa/db.index msa/db.ffindex
        elif [ -e msa/db && -e msa/db.ffdata ]
        then
            true
        else
            echo "An msa was specified but there was no db.index or db.ffindex file!" 1>&2
            exit 1
        fi


        # split the db
        TARGET_CLUSTER_SIZE=10000
        NCLUSTERS=\$(wc -l < "msa/db.index")
        NSPLITS=\$(( (\${NCLUSTERS} + \${TARGET_CLUSTER_SIZE} + 1) / \${TARGET_CLUSTER_SIZE} ))

        mmseqs splitdb "msa/db" "tmp_split_msa" --split \${NSPLITS}

        for f in split_msa_*.index
        do
            BASENAME="${f%.index}"
            DIRNAME="${BASENAME#tmp_}"

            mkdir "${DIRNAME}"
            mv "${f}" "${DIRNAME}/db.index"
            ln -s "${DIRNAME}/db.index" "${DIRNAME}/db.ffindex"

            mv "${BASENAME}" "${DIRNAME}/db"
            ln -s "${DIRNAME}/db" "${DIRNAME}/db.ffdata"
        done
        """

    }

} else {

    msa = combinedMSAs
    splitMSAs4EnrichMSA = splitMSAs4DecideIfUser

}


//
// STEP 4 - Get a tree for each cluster MSA
//


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
    run_remote

    input:
    file "msa" from splitMSAs4EnrichMSA
    file "enrich" from enrichdb

    output:
    file "fas" into splitFas

    script:
    """
    mkdir profile
    mmseqs msa2profile msa/db profile/db --match-mode 1 --match-ratio 1

    mkdir search tmp
    mmseqs search profile/db enrich/db search/db tmp -a

    mkdir fas
    mmseqs result2msa profile/db enrich/db search/db fas/db

    mv fas/db fas/db.ffdata
    mv fas/db.dbtype fas/db.ffdata.dbtype
    mv fas/db.index fas/db.ffindex

    rm -rf -- profile search tmp
    """
}


/*
 * Create HHsuite databases for each MSA chunk.
 */
process fasToHHDB {
    label "hhsuite"
    label "big_task"

    when:
    run_remote

    input:
    file "fas" from splitFas
    file "hhdata" from hhdata

    output:
    file "hhdb" into splitHHDB

    script:
    """
    mkdir -p hhdb
    cp fas/db.ffdata hhdb/db_fasta.ffdata
    cp fas/db.ffdata.dbtype hhdb/db_fasta.ffdata.dbtype
    cp fas/db.ffindex hhdb/db_fasta.ffindex

    mpirun -np "${task.cpus}" ffindex_apply_mpi \
        hhdb/db_fasta.ff{data,index} \
        -i hhdb/db_a3m.ffindex \
        -d hhdb/db_a3m.ffdata \
        -- \
        run_fas_to_a3m.sh

    mpirun -np "${task.cpus}" cstranslate_mpi \
        -i hhdb/db_a3m \
        -o hhdb/db_cs219 \
        -x 0.3 \
        -c 4 \
        -b \
        -I a3m \
        -A hhdata/cs219.lib \
        -D hhdata/context_data.lib

    mpirun -np "${task.cpus}" ffindex_apply_mpi \
        hhdb/db_a3m.ff{data,index} \
        -i hhdb/db_hhm.ffindex \
        -d hhdb/db_hhm.ffdata \
        -- \
        hhmake -i stdin -o stdout -v 0
    """
}


splitHHDB.into {
    splitHHDB4CombineSplitHHDBs;
    splitHHDB4SearchSelf;
    splitHHDB4SearchPfam;
    splitHHDB4SearchScop;
    splitHHDB4SearchPdb;
    splitHHDB4SearchUniref;
}


/*
 * Combine the HHsuite database chunks into a final database.
 * We can use this database to search against.
 */
process combineSplitHHDBs {

    label "ffdb"
    label "small_task"

    publishDir "${params.outdir}"

    when:
    run_remote

    input:
    file "split_hhdata_*" from splitHHDB4CombineSplitHHDBs.collect()

    output:
    file "hhdata" into HHDB

    script:
    """
    mkdir unsorted hhdata

    ffdb combine \
      -d unsorted/db_fasta.ffdata \
      -i unsorted/db_fasta.ffindex \
      split_hhdata_*/db_fasta.ff{data,index}

    ffdb combine \
      -d unsorted/db_a3m.ffdata \
      -i unsorted/db_a3m.ffindex \
      split_hhdata_*/db_a3m.ff{data,index}

    ffdb combine \
      -d unsorted/db_cs219.ffdata \
      -i unsorted/db_cs219.ffindex \
      split_hhdata_*/db_cs219.ff{data,index}

    ffdb combine \
      -d unsorted/db_hhm.ffdata \
      -i unsorted/db_hhm.ffindex \
      split_hhdata_*/db_hhm.ff{data,index}

    sort -k3 -n unsorted/db_cs219.ffindex | cut -f1 > sorting.dat

    ffindex_order sorting.dat unsorted/db_cs219.ff{data,index} hhdata/db_cs219.ff{data,index}
    ffindex_order sorting.dat unsorted/db_a3m.ff{data,index} hhdata/db_a3m.ff{data,index}
    ffindex_order sorting.dat unsorted/db_fasta.ff{data,index} hhdata/db_fasta.ff{data,index}
    ffindex_order sorting.dat unsorted/db_hhm.ff{data,index} hhdata/db_hhm.ff{data,index}

    # rm -rf -- unsorted
    """
}


//
// STEP
//


process searchSelf {

    label "hhsuite"
    label "big_task"

    when:
    run_remote && !params.hhmatches_self

    input:
    file "subdb" from splitHHDB4SearchSelf
    file "hhdb" from HHDB

    output:
    file "results" into splitSelfMatchesTmp

    script:
    """
    mkdir results
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i subdb/db_hhm \
      -d hhdb/db \
      -o results/db_hhr \
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


process searchPfam {

    label "hhsuite"
    label "big_task"

    when:
    run_remote && params.hhmatches_pfam

    input:
    file "hhdata" from splitHHDB4SearchPfam
    file "pfam" from pfamdb

    output:
    file "results" into splitPfamMatchesTmp

    script:
    """
    mkdir results
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i hhdata/db_hhm \
      -d pfam/pfam \
      -o results/db_hhr \
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

    splitPfamMatches = splitSelfMatchesTmp

}


process searchScop {

    label "hhsuite"
    label "big_task"

    when:
    run_remote && !params.hhmatches_scop

    input:
    file "hhdata" from splitHHDB4SearchScop
    file "scop" from scopdb

    output:
    file "results" into splitScopMatchesTmp

    script:
    """
    mkdir results
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i hhdata/db_hhm \
      -d scop/scop \
      -o results/db \
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


process searchPdb {

    label "hhsuite"
    label "big_task"

    when:
    run_remote && !params.hhmatches_pdb

    input:
    file "hhdata" from splitHHDB4SearchPdb
    file "pdb" from pdbdb

    output:
    file "results" into splitPdbMatchesTmp

    script:
    """
    mkdir results
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i hhdata/db_hhm \
      -d pdb/pdb \
      -o results/db \
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


process searchUniref {

    label "hhsuite"
    label "big_task"

    when:
    run_remote && !params.hhmatches_uniref

    input:
    file "hhdata" from splitHHDB4SearchUniref
    file "uniref" from unirefdb

    output:
    file "results" into splitUnirefMatchesTmp

    script:
    """
    mkdir results
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i hhdata/db_hhm \
      -d uniref/uniref \
      -o results/db \
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

    splitUnirefMatches = splitSelfMatchesTmp

}


splitSelfMatches.collect().map { ds -> ["self", ds] }
    .mix(
        splitPfamMatches.collect().map { ds -> ["pfam", ds] },
        splitScopMatches.collect().map { ds -> ["scop", ds] },
        splitPdbMatches.collect().map { ds -> ["pdb", ds] },
        splitUnirefMatches.collect().map { ds -> ["uniref", ds] },
    )
    .set { splitMatches }


process combineSearchMatches {

    label "ffdb"
    label "small_task"

    publishDir "${params.outdir}/hhmatches"

    tag "${database}"

    input:
    set val(database), file("split_results_*") from splitMatches

    output:
    set val(database), file("${database}_matches") into combinedMatches

    script:
    """
    mkdir -p "${database}"

    ffdb combine \
      -d "${database}/db_hhr.ffdata" \
      -i "${database}/db_hhr.ffindex" \
      split_results_*/db_hhr.ff{data,index}
    """
}

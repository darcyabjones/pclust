#!/usr/bin/env nextflow


//
// Database creation
//

process create_db {

    label 'mmseqs'
    label "small_task"

    input:
    path "seqs.fasta"

    output:
    path "seqdb", emit: db

    script:
    """
    mkdir -p "seqdb"
    mmseqs createdb "seqs.fasta" "seqdb/db" --max-seq-len 14000
    """
}


process msa_to_profile {

    label "mmseqs"
    label "small_task"

    input:
    path "msas"

    output:
    path "profiles", emit: profiles

    script:
    """
    mkdir -p "profiles"

    mmseqs msa2profile \
      "msas/db" \
      "profile/db" \
      --match-mode 1 \
      --match-ratio 1
    """
}


process create_profiles_from_search {

    label "mmseqs"
    label "big_task"

    input:
    path "query_db"
    path "target_db"
    path "matches"

    output:
    path "profiles", emit: profiles

    script:
    """
    mkdir -p "profiles"

    mmseqs result2profile \
      "query_db/db" \
      "target_db/db" \
      "matches/db" \
      "profiles/db" \
      --threads "${task.cpus}"
    """
}


//
// Searching
//


/*
 * Quickly search a database.
 * No alignment information, just matches.
 */
process search_profiles_quick {

    label "mmseqs"
    label "big_task"

    input:
    path "profiles"
    path "db"

    output:
    path "matches", emit: matches

    script:
    """
    mkdir -p "tmp"
    mkdir -p "matches"

    mmseqs search \
      "profiles/db" \
      "db/db" \
      "matches/db" \
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
 * Uses higher -s parameter and computes alignment info.
 */
process search_profiles_sensitive {

    label "mmseqs"
    label "big_task"

    input:
    path "profiles"
    path "db"

    output:
    path "matches", emit: matches

    script:
    """
    mkdir -p "matches" "tmp"

    mmseqs search \
      "profiles/db" \
      "db/db" \
      "matches/db" \
      "tmp" \
      --threads "${task.cpus}" \
      -s 6 \
      -a \
      -e 0.00001 \
      --db-load-mode 0 \
      --split 0 \
      --split-mode 1

    rm -rf -- tmp
    """
}


/*
 * Search the cluster profiles against the cluster consensus sequences.
 */
process search_profiles_against_self {

    label 'mmseqs'
    label "big_task"

    input:
    path "profiles"

    output:
    path "matches", emit: matches

    script:
    """
    mkdir "tmp"
    mkdir "matches"

    mmseqs search \
      "profiles/db" \
      "profiles/db_consensus" \
      "matches/db" \
      "tmp" \
      --threads "${task.cpus}" \
      --max-seqs 300 \
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


//
// Enrichment
//


process enrich_msas_from_search {

    label "mmseqs"
    label "small_task"

    input:
    path "msas"
    path "profiles"
    path "seqs"
    path "matches"

    output:
    path "enriched", emit: enriched

    script:
    """
    mkdir "match_msas"
    mmseqs result2msa \
      "profiles/db" \
      "seqs/db" \
      "matches/db" \
      "match_msas/db" \
      --skip-query 1

    mkdir "concatenated"
    ffdb join_concat \
      -d "concatenated/db" \
      -i "concatenated/db.index" \
      "msas/db" "match_msas/db" \
      "msas/db.index" "match_msas/db.index"

    mkdir "enriched"
    mpirun -np ${task.cpus} mmseqs apply \
        "concatenated/db" \
        "enriched/db" \
        --threads 1 \
        -- \
        tidy_joined_alignments.py

    # Would be good to sort by num msa columns to help balance load.
    # Not necessary short-term, order seems like that anyway.

    cp -L "match_msas/db.dbtype" "enriched/db.dbtype"

    ORIG="\${PWD}"
    cd "enriched"
    ln -s "db" "db.ffdata"
    ln -s "db.dbtype" "db.ffdata.dbtype"
    ln -s "db.index" "db.ffindex"

    cd "\${ORIG}"

    rm -rf -- "concatenated" "match_msas"
    """
}



//
// Clustering
//

/*
 * Perform basic mmseqs workflow.
 */
process cascade_cluster {

    label 'mmseqs'
    label "big_task"

    input:
    path "seqs"

    output:
    path "clusters", emit: clusters

    script:
    """
    mkdir -p "clusters"
    mkdir -p "tmp"

    mmseqs cluster \
      "seqs/db" \
      "clusters/db" \
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
process extract_cluster_stats {

    label 'mmseqs'
    label "big_task"

    input:
    path "seqs"
    path "clusters"

    output:
    path "clusters.tsv", emit: clusters
    path "clusters_rep.fasta", emit: representative
    path "clusters_stats.tsv", emit: stats

    script:
    """
    mmseqs createtsv \
      "seqs/db" \
      "seqs/db" \
      "clusters/db" \
      "clusters.tsv" \
      --threads ${task.cpus}

    sed -i '1i cluster\tmember' "clusters.tsv"


    mmseqs result2repseq \
      "seqs/db" \
      "clusters/db" \
      "clusters_rep" \
      --threads ${task.cpus}

    mmseqs result2flat \
      "seqs/db" \
      "seqs/db" \
      "clusters_rep" \
      "clusters_rep.fasta" \
      --use-fasta-header


    # Do an all vs centroid alignment for each cluster.
    mmseqs align \
      "seqs/db" \
      "seqs/db" \
      "clusters/db" \
      align \
      -a \
      --threads ${task.cpus}

    mmseqs convertalis \
      "seqs/db" \
      "seqs/db" \
      align \
      "clusters_stats.tsv" \
      --threads ${task.cpus} \
      --format-mode 0 \
      --format-output 'query,target,evalue,qcov,tcov,gapopen,pident,nident,mismatch,raw,bits,qstart,qend,tstart,tend,qlen,tlen,alnlen'

    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' "clusters_stats.tsv"
    """
}


/*
 * Cluster the cluster profiles based on the profile/consensus search.
 */
process cluster_profile_self_matches {

    label 'mmseqs'
    label "big_task"

    input:
    path "profiles"
    path "matches"

    output:
    path "clusters", emit: clusters

    script:
    """
    mkdir -p "clusters"

    mmseqs clust \
      "profiles/db" \
      "matches/db" \
      "clusters/db" \
      --threads "${task.cpus}" \
      --cluster-mode 2
    """
}


/*
 * Merge the cascade and profile clustering results to get final clusters.
 */
process merge_clusters_by_other_clusters {

    label "mmseqs"
    label "small_task"

    input:
    path "seqs"
    path "child_clusters"
    path "parent_clusters"

    output:
    path "clusters", emit: clusters

    script:
    """
    mkdir "clusters"

    mmseqs mergeclusters \
      "seqs/db" \
      "clusters/db" \
      "child_clusters/db" \
      "parent_clusters/db"
    """
}


process get_cluster_seqs {

    label "mmseqs"
    label "small_task"

    input:
    path "seqs"
    path "clusters"

    output:
    path "cluster_seqs", emit: "seqs"

    script:
    """
    mkdir -p "cluster_seqs"

    mmseqs createseqfiledb \
      "seqs/db" \
      "clusters/db" \
      "cluster_seqs/db"

    ORIG="\${PWD}"
    cd cluster_seqs
    ln -s db db.ffdata
    ln -s db.dbtype db.ffdata.dbtype
    ln -s db.index db.ffindex
    cd \${ORIG}
    """
}


//
// Workflows
//

workflow enrich_profiles {

    get:
    profiles
    enrichdb

    main:
    matches = search_profiles_quick(profiles, enrichdb)
    enriched = create_profiles_from_search(profiles, enrichdb, matches)

    emit:
    enriched
}


workflow profile_cluster {

    get:
    seqs
    clusters
    enrich // Should be a bool value channel
    enrichdb

    main:
    profiles = create_profiles_from_search(seqs, seqs, clusters)

    if ( enrich ) {
        profiles_to_use = enrich_profiles(profiles, enrichdb)
    } else {
        profiles_to_use = profiles
    }

    self_matches = search_profiles_against_self(profiles_to_use)
    self_matches_clusters = cluster_profile_self_matches(profiles_to_use, self_matches)
    profile_clusters = merge_clusters_by_other_clusters(seqs, clusters, self_matches_clusters)
    (pc_tsv, pc_representative, pc_stats) = extract_cluster_stats(seqs, profile_clusters)

    emit:
    profile_clusters
    pc_tsv
    pc_representative
    pc_stats
}


workflow cluster {

    get:
    seqs
    enrich // Should be a bool value channel
    enrichdb // Must exist if `enrich` is true

    main:
    cc = cascade_cluster(seqs)
    (cc_tsv, cc_representative, cc_stats) = extract_cluster_stats(seqs, cc)
    (pc, pc_tsv, pc_representative, pc_stats) = profile_cluster(seqs, cc, enrichdb, enrich)

    // pc_seqs = get_cluster_seqs(seqs, pc)

    emit:
    cc
    cc_tsv
    cc_representative
    cc_stats
    pc
    pc_tsv
    pc_representative
    pc_stats
}


/*
 * Enrich the MSAs by searching against a database.
 * This will give extra information that might not be there
 * because we have quite strict coverage requirements.
 * It will help matching domains later on.
 *
 * We use MMSeqs because it is MUCH faster than the HHBlits enrichment
 * method.
 */
workflow enrich_msas {

    get:
    msas
    enrichdb

    main:
    profiles = msa_to_profile(msas)
    profile_matches = search_profiles_sensitive(profiles, enrichdb)
    enriched = enrich_msas_from_search(msas, profiles, enrichdb, profile_matches)

    emit:
    enriched
}

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
    path "msa_profiles", emit: profiles

    script:
    """
    mkdir -p "msa_profiles"

    cp -rL msas msas_tmp
    awk 'BEGIN { printf("%c%c%c%c",11,0,0,0); exit; }' > msas_tmp/db.dbtype

    mmseqs msa2profile \
      "msas_tmp/db" \
      "msa_profiles/db" \
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
    path "search_profiles", emit: profiles

    script:
    """
    mkdir -p "search_profiles"

    mmseqs result2profile \
      "query_db/db" \
      "target_db/db" \
      "matches/db" \
      "search_profiles/db" \
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
      -c 0.7 \
      --cov-mode 0 \
      --min-seq-id 0.1 \
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
    path "enriched_msas", emit: enriched

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

    mkdir "enriched_msas"
    mpirun -np ${task.cpus} mmseqs apply \
        "concatenated/db" \
        "enriched_msas/db" \
        --threads 1 \
        -- \
        tidy_joined_alignments.py

    # Would be good to sort by num msa columns to help balance load.
    # Not necessary short-term, order seems like that anyway.

    cp -L "match_msas/db.dbtype" "enriched_msas/db.dbtype"

    ORIG="\${PWD}"
    cd "enriched_msas"
    ln -s "db" "db.ffdata"
    ln -s "db.dbtype" "db.ffdata.dbtype"
    ln -s "db.index" "db.ffindex"

    cd "\${ORIG}"

    # rm -rf -- "concatenated" "match_msas"
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
    path "cascade", emit: clusters

    script:
    """
    mkdir -p "cascade"
    mkdir -p "tmp"

    mmseqs cluster \
      "seqs/db" \
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
process extract_cluster_stats {

    label 'mmseqs'
    label "big_task"

    input:
    path "seqs"
    path "clusters"
    val name

    output:
    path "${name}.tsv", emit: clusters
    path "${name}_rep.fasta", emit: representative
    path "${name}_stats.tsv", emit: stats

    script:
    """
    mmseqs createtsv \
      "seqs/db" \
      "seqs/db" \
      "clusters/db" \
      "${name}.tsv" \
      --threads ${task.cpus}

    sed -i '1i cluster\tmember' "${name}.tsv"


    mmseqs result2repseq \
      "seqs/db" \
      "clusters/db" \
      "clusters_rep" \
      --threads ${task.cpus}

    mmseqs result2flat \
      "seqs/db" \
      "seqs/db" \
      "clusters_rep" \
      "${name}_rep.fasta" \
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
      "${name}_stats.tsv" \
      --threads ${task.cpus} \
      --format-mode 0 \
      --format-output 'query,target,evalue,qcov,tcov,gapopen,pident,nident,mismatch,raw,bits,qstart,qend,tstart,tend,qlen,tlen,alnlen'

    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' "${name}_stats.tsv"
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
    path "profile", emit: clusters

    script:
    """
    mkdir "profile"

    mmseqs mergeclusters \
      "seqs/db" \
      "profile/db" \
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
    val name

    output:
    path "${name}_seqs", emit: "seqs"

    script:
    """
    mkdir -p "${name}_seqs"

    mmseqs createseqfiledb \
      "seqs/db" \
      "clusters/db" \
      "${name}_seqs/db"

    # The mmseqs createseqfiledb type is Generic, which doesn't play well
    # with splitdb. This makes the type Protein.
    cp seqs/db.dbtype "${name}_seqs/db.dbtype"
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
    (pc_tsv, pc_representative, pc_stats) = extract_cluster_stats(seqs, profile_clusters, "profile")

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
    (cc_tsv, cc_representative, cc_stats) = extract_cluster_stats(seqs, cc, "cascade")
    (pc, pc_tsv, pc_representative, pc_stats) = profile_cluster(seqs, cc, enrich, enrichdb)

    //pc_seqs = get_cluster_seqs(seqs, pc, "profile")

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

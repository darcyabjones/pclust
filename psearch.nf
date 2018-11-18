#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    =================================
    pclust/pupdate
    =================================

    Usage:

    abaaab

    Mandatory Arguments:
      --proteins              description
      --proteins_db
      --global_profile
      --global_clusters
      --global_seqs

    Options:
      --trees
      --nomsa

    Outputs:

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

params.query_profile = false
params.query_clusters = false
params.query_seqs = false
params.target_seqs = false
params.target_db = false



if ( !params.query_profile || !(params.query_clusters && params.query_seqs) ) {
    log.info "Need either the query profile or query clusters + seqdb"
    exit 1
}


if (params.target_db) {
    targetSeqs = Channel.fromPath( params.target_db )
} else if (params.target_seqs) {
    proteins = Channel.fromPath( params.target_seqs )
    
    /*
     * Create the mmseqs2 sequence database
     */
    process createSequenceDB {
        label 'mmseqs'
        publishDir "${params.outdir}/sequences"
    
        input:
        file fasta from proteins
    
        output:
        file "proteins" into targetSeqs
    
        """
        mkdir -p proteins
        mmseqs createdb "${fasta}" proteins/db --max-seq-len 14000
        """
    }
} else {
    log.info "Need either --target_seqs or --target_db"
    exit 1
}


if ( params.query_profile ) {
    queryProfile = Channel.fromPath( params.query_profile )
} else {
    queryClusters = Channel.fromPath( params.query_clusters )
    querySeqs = Channel.fromPath( params.query_seqs )

    /*
     * Create profile to search against target.
     */
    process createGlobalProfile {
        label "mmseqs"
        publishDir "${params.outdir}/search"
    
        input:
        file "clusters" from queryClusters
        file "seqs" from querySeqs
    
        output:
        file "global_profile" into queryProfile
    
        """
        mkdir -p global_profile
    
        # Create profiles for each cluster.
        # Generates the profile_consensus file too.
        mmseqs result2profile \
          seqs/db \
          seqs/db \
          clusters/db \
          global_profile/db \
          --threads ${task.cpus}
        """
    }
}


/*
 * Search the query profile against the target seqs.
 */
process searchTarget {
    label "mmseqs"
    publishDir "${params.outdir}/search"

    input:
    file "profile" from queryProfile
    file "seqs" from targetSeqs

    output:
    file "search" into searchResults 
    file "search.tsv" into searchResultsTsv

    """
    mkdir -p tmp
    mkdir -p search
    mmseqs search \
      profile/db \
      seqs/db \
      search/db \
      tmp \
      -a \
      -s 6.0 \
      --max-seqs 5000 \
      --e-profile 0.01 \
      --rescore-mode 1 \
      --num-iterations 3

    mmseqs convertalis \
      profile/db \
      seqs/db \
      search/db \
      "search.tsv" \
      --threads ${task.cpus} \
      --format-mode 0 \
      --format-output "query target evalue qcov tcov gapopen pident nident mismatch raw bits qstart qend tstart tend qlen tlen alnlen cigar"

    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen\tcigar' search.tsv
    """
}


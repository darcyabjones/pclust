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

params.global_profile = false
params.global_clusters = false
params.global_seqs = false
params.proteins = false
params.proteins_db = false



if (params.global_profile) {
    globalProfile = Channel.fromPath( params.global_profile )
} else if ( params.global_clusters && params.global_seqs ) {
    globalClusters = Channel.fromPath( params.global_clusters )
    globalSeqs = Channel.fromPath( params.global_seqs )
} else {
    log.info "Need either the global profile or global clusters + seqdb"
    exit 1
}

if (params.proteins_db) {
    seq = Channel.fromPath( params.proteins_db )
} else if (params.proteins) {
    proteins = Channel.fromPath( params.proteins )
    
    /*
     * Create the mmseqs2 sequence database
     */
    process createSequenceDB {
        label 'mmseqs'
        publishDir "sequences"
    
        input:
        file fasta from proteins
    
        output:
        file "proteins" into seq
    
        """
        mkdir -p proteins
        mmseqs createdb "${fasta}" proteins/db --max-seq-len 14000
        """
    }
} else {
    log.info "Need either --proteins or --proteins_db"
    exit 1
}


if ( !params.global_profile ) {
    process createGlobalProfile {
        label "mmseqs"
        publishDir "clusters"
    
        input:
        file "clusters" from globalClusters
        file "seqs" from globalSeqs
    
        output:
        file "global_profile" into globalProfile
    
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


process searchGlobal {
    label "mmseqs"
    publishDir "clusters"

    input:
    file "seqs" from seq
    file "profile" from globalProfile

    output:
    file "search" into globalSearchResults 
    file "search.tsv" into globalSearchResultsTsv

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



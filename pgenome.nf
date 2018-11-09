#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    =================================
    pclust/pgenome
    =================================

    Usage:

    abaaab

    Mandatory Arguments:
      --genomes
      --gffs
      --clusters
      --proteindb               description

    Options:
      --non-existant          description

    Outputs:

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

params.gffs = false
if ( params.gffs ) {
    gffs = Channel
        .fromPath( params.gffs )
        .map { [it.baseName, it] }
}

params.genomes = Channel
    .fromPath( params.genomes )
    .map { [it.baseName, it] }

clusters = Channel.fromPath( params.clusters )
seq = Channel.fromPath( params.proteindb )


/*
 * Construct genomes database
 */
process createGenomeSequenceDB {
    label 'mmseqs'
    publishDir "genomes"
    tag { label }

    input:
    set val(label), file(fasta) from combinedGenomeFasta

    output:
    set val(label), file("${label}") into genomeSeq

    """
    mkdir -p "${label}"
    mmseqs createdb "${fasta}" "${label}/db" --dont-split-seq-by-len
    """
}


/*
 * Search for the clusters in the original genome sequences.
 */
process searchGenomes {
    label 'mmseqs'
    publishDir "genomes"
    tag { label }

    input:
    set val(label), file("genome") from genomeSeq
    file "clusters" from clusters
    file "seqs" from seq

    output:
    set val(label), file("${label}.tsv") into searchResults

    """
    mkdir -p tmp

    # Create profiles for each cluster.
    mmseqs result2profile \
      seqs/db \
      seqs/db \
      clusters/db \
      profile \
      --threads ${task.cpus}

    # Search profile vs genome
    # Search parameters are slightly more conservative than default.
    mmseqs search \
      profile \
      genome/db \
      result \
      tmp \
      --threads ${task.cpus} \
      --realign \
      --gap-open 15 \
      --gap-extend 2 \
      --cov-mode 2 \
      --rescore-mode 2 \
      --min-length 20 \
      --orf-start-mode 1 \
      --use-all-table-starts true

    # Extract matches from results
    mmseqs convertalis \
      profile \
      genome/db \
      result "${label}.tsv" \
      --threads ${task.cpus} \
      --format-mode 0 \
      --format-output "query target evalue qcov tcov gapopen pident nident mismatch raw bits qstart qend tstart tend qlen tlen alnlen cigar qframe tframe"

    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen\tcigar\tqframe\ttframe' profile_matches.tsv
    """
}

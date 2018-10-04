#!/usr/bin/env nextflow

params.fastas = "$baseDir/data/*.faa"

fastas = Channel.fromPath( params.fastas )


// Combine the fasta files in preparation for clustering
// Need step to rename fastas to include filename.
process combineFasta {

    input:
    file "*.fa*" from fastas.collect()

    output:
    file "combined.faa" into combinedFasta

    """
    cat *.fa* > combined.faa
    """
}


process createSequenceDB {
    container "soedinglab/mmseqs2"

    input:
    file fasta from combinedFasta

    output:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  into sequenceDB

    """
    mmseqs createdb "${fasta}" sequence --max-seq-len 14000
    """
}


sequenceDB.into { sequenceDB1; sequenceDB2; sequenceDB3; sequenceDB4; sequenceDB5; sequenceDB6; sequenceDB7; sequenceDB8; sequenceDB9; sequenceDB10; sequenceDB11 }

process mmseqsClust {
    container "soedinglab/mmseqs2"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB1

    output:
    set "clusters", "clusters.index" into clustDB

    """
    mkdir -p tmp
    mmseqs cluster sequence clusters tmp -s 5 --cluster-steps 5 --min-seq-id 0.3
    """
}

clustDB.into {clustDB1; clustDB2; clustDB3 }

process mmseqsExtractClusters {
    container "soedinglab/mmseqs2"

    publishDir "clusters"

    input:
    set "clusters", "clusters.index" from clustDB1
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB2

    output:
    file "clusters.tsv" into uniclust30TSV
    file "clusters_rep.fasta" into uniclust30Fasta
    file "align.m8" into align

    """
    mmseqs createtsv \
      sequence \
      sequence \
      clusters \
      clusters.tsv

    mmseqs result2repseq \
      sequence \
      clusters \
      clusters_rep

    mmseqs result2flat \
      sequence \
      sequence \
      clusters_rep \
      clusters_rep.fasta \
      --use-fasta-header

    mmseqs align sequence sequence clusters align -a
    mmseqs convertalis sequence sequence align align.m8
    """
}


process mmseqsProfileClust {
    container "soedinglab/mmseqs2"

    input:
    set "clusters", "clusters.index" from clustDB2
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB3

    output:
    set "profile_clusters", "profile_clusters.index" into profileClustDB

    """
    mmseqs result2profile sequence sequence clusters profile

    mkdir -p tmp
    mmseqs search profile profile_consensus result tmp -s 7.5 -e 0.05 --add-self-matches --num-iterations 3
    mmseqs clust profile result profile_clusters_consensus

    mmseqs mergeclusters sequence profile_clusters clusters profile_clusters_consensus
    """
}


process mmseqsExtractProfileClusters {
    container "soedinglab/mmseqs2"

    publishDir "profile_clusters"

    input:
    set "profile_clusters", "profile_clusters.index" from profileClustDB
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB4

    output:
    file "clusters.tsv" into uniclust30TSV2
    file "clusters_rep.fasta" into uniclust30Fasta2
    file "align.m8" into align2

    """
    mmseqs createtsv \
      sequence \
      sequence \
      profile_clusters \
      clusters.tsv

    mmseqs result2repseq \
      sequence \
      profile_clusters \
      clusters_rep

    mmseqs result2flat \
      sequence \
      sequence \
      clusters_rep \
      clusters_rep.fasta \
      --use-fasta-header

    mmseqs align sequence sequence profile_clusters align -a
    mmseqs convertalis sequence sequence align align.m8
}


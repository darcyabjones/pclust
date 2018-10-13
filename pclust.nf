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
      --genomes               description

    Options:
      --non-existant          description

    Outputs:

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

params.genomes = "$baseDir/data/*.{fasta,gff3}"
genomes = Channel.fromFilePairs( params.genomes, flat: true )
genomes.tap { genomes1 }


/*
 * Extract protein sequences from genomes and GFF.
 */
process extractProteins {
    container "quay.io/biocontainers/genometools-genometools:1.5.10--h470a237_1"

    input:
    set val(label), file(fasta), file(gff) from genomes1

    output:
    file "${label}.faa" into proteins

    """
    gt extractfeat \
      -type CDS \
      -translate \
      -matchdescstart \
      -join \
      -retainids \
      -seqfile "${fasta}" \
      "${gff}" \
      > "${label}.faa"
    """
}

proteins.into { proteins1; proteins2 }


/*
 * Combine all proteins into a single fasta file.
 */
process combineFasta {

    input:
    file "*.fa*" from proteins1.collect()

    output:
    file "combined.faa" into combinedFasta

    """
    cat *.fa* > combined.faa
    """
}


/*
 * Create the mmseqs2 sequence database
 */
process createSequenceDB {
    container "soedinglab/mmseqs2"
    publishDir "clusters"
    label 'mmseqs'

    input:
    file fasta from combinedFasta

    output:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  into sequenceDB

    """
    mmseqs createdb "${fasta}" sequence --max-seq-len 14000
    """
}

sequenceDB.into {
    sequenceDB1; 
    sequenceDB2;
    sequenceDB3;
    sequenceDB4;
    sequenceDB5;
    sequenceDB6;
    sequenceDB7
    }


/*
 * Do the first clustering pass.
 * This tends to get clusters down to approximately 50% identity.
 */
process clust {
    container "soedinglab/mmseqs2"
    publishDir "clusters"
    label 'mmseqs'

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
         "sequence_h", "sequence_h.index"  from sequenceDB1

    output:
    set "clusters", "clusters.index" into clustDB

    """
    mkdir -p tmp
    mmseqs cluster sequence clusters tmp -s 5 --cluster-steps 5 --min-seq-id 0.3
    """
}

clustDB.into {
    clustDB1;
    clustDB2;
    clustDB3;
    clustDB4;
    clustDB5
}


/*
 * Perform the second clustering pass using sequence profiles.
 * This gets clusters down to about 10%-20% identity.
 */
process profileClust {
    container "soedinglab/mmseqs2"
    publishDir "clusters"
    label 'mmseqs'

    input:
    set "clusters", "clusters.index" from clustDB1
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
        "sequence_h", "sequence_h.index"  from sequenceDB2

    output:
    set "profile_clusters", "profile_clusters.index" into profileClustDB

    """
    # Create profiles for each cluster.
    # Generates the profile_consensus file too.
    mmseqs result2profile sequence sequence clusters profile

    # Search the profiles against the profile consensus sequences.
    # Uses an iterative strategy.
    mkdir -p tmp
    mmseqs search \
      profile \
      profile_consensus \
      result \
      tmp \
      -s 7.5 \
      -e 0.05 \
      --add-self-matches \
      --num-iterations 3

    # Cluster the matches of profiles vs consensus sequences.
    mmseqs clust profile result profile_clusters_consensus

    # Merge the original clusters with the new ones to get all sequences back.
    mmseqs mergeclusters \
      sequence \
      profile_clusters \
      clusters \
      profile_clusters_consensus
    """
}

profileClustDB.into { profileClustDB1; profileClustDB2 }


/*
 * Merge the two clustering results into a single channel that we can look at.
 */
allClusters = Channel
    .value("identity")
    .combine(clustDB2)
    .concat(
        Channel.value("profile").combine( profileClustDB2 )
    )

allClusters.into { allClusters1; allClusters2; allClusters3 }


/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractClusterStats {
    container "soedinglab/mmseqs2"
    publishDir { "clusters/${type}" }
    tag { type }
    label 'mmseqs'

    input:
    set val(type), file(clusters), file(clusters_index) from allClusters1
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
        "sequence_h", "sequence_h.index" from sequenceDB4

    output:
    file "${clusters}.tsv" into clustersTSV
    file "${clusters}_rep.fasta" into clustersRepFasta
    file "${clusters}_stats.tsv" into clustersStats

    """
    mmseqs createtsv \
      sequence \
      sequence \
      "${clusters}" \
      "${clusters}.tsv"

    mmseqs result2repseq \
      sequence \
      "${clusters}" \
      "${clusters}_rep"

    mmseqs result2flat \
      sequence \
      sequence \
      "${clusters}_rep" \
      "${clusters}_rep.fasta" \
      --use-fasta-header

    # Do an all vs centroid alignment for each cluster.
    mmseqs align sequence sequence "${clusters}" align -a
    mmseqs convertalis sequence sequence align "${clusters}_stats.tsv"
    """
}


/*
 * Extract sequences from clusters.
 */
process extractClusterFastaLikes {
    container "soedinglab/mmseqs2"
    tag { type }
    label 'mmseqs'

    input:
    set val(type), file(clusters), file(clusters_index) from allClusters2
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
        "sequence_h", "sequence_h.index" from sequenceDB5

    output:
    set val(type), file("${clusters}_sequences.fastalike") into allClustersFastaLikes

    """
    mmseqs createseqfiledb sequence "${clusters}" "${clusters}_sequences"
    mmseqs result2flat \
        sequence \
        sequence \
        "${clusters}_sequences" \
        "${clusters}_sequences.fastalike"
    """
}


/*
 * Extract the infividual fasta sequences from the fasta-like file.
 */
process extractClusterFastas {

    tag { type }

    input:
    set val(type), file(fastalike) from allClustersFastaLikes

    output:
    set val(type), file("*.fasta") into allClustersFasta

    """
    extract_fastalike.py "${fastalike}"
    """
}

allClustersFasta.transpose().into{ allClustersFasta1; allClustersFasta2 }

process mafft {
    container "pclust/mafft_mmseqs2"
    publishDir { "msas/mafft/${type}"}
    tag { type }

    input:
    set val(type), file(fasta) from allClustersFasta1

    output:
    set val(type), file("${fasta.baseName}.faa") into mafftMsas

    """
    mafft \
        --amino \
        --thread 2 \
        --retree 2 \
        --maxiterate 1 \
        "${fasta}" \
    > "${fasta.baseName}.faa"
    """
}

process muscle {
    container "quay.io/biocontainers/muscle:3.8.1551--h2d50403_3"
    publishDir { "msas/muscle/${type}"}
    tag { type }

    input:
    set val(type), file(fasta) from allClustersFasta2

    output:
    set val(type), file("${fasta.baseName}.faa") into muscleMsas

    """
    muscle \
      -in "${fasta}" \
      -out "${fasta.baseName}.faa" \
      -maxiters 2 \
      -seqtype protein
    """
}

/*
 * Extract sequences from clusters.
 */
process mmseqs2MSA {
    container "soedinglab/mmseqs2"
    publishDir { "msas/mmseqs/${type}"}
    tag { type }
    label 'mmseqs'

    input:
    set val(type), file(clusters), file(clusters_index) from allClusters3
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
        "sequence_h", "sequence_h.index" from sequenceDB3

    output:
    set val(type), file("msa.fastalike") into msaFastaLikes

    """
    mmseqs result2msa sequence sequence "${clusters}" msa
    mmseqs result2flat sequence sequence msa msa.fastalike
    """
}

/*
 * Extract the infividual fasta sequences from the fasta-like file.
 */
process extractMSAFastas {
    publishDir { "msas/mmseqs/${type}"}
    tag { type }

    input:
    set val(type), file(fastalike) from msaFastaLikes

    output:
    set val(type), file("*.fasta") into mmseqsMsas

    """
    extract_fastalike.py "${fastalike}"
    """
}

process muscleRefine {
    container "quay.io/biocontainers/muscle:3.8.1551--h2d50403_3"
    publishDir { "msas/muscle_refine/${type}"}
    tag { type }

    input:
    set val(type), file(fasta) from mmseqsMsas.transpose()

    output:
    set val(type), file("${fasta.baseName}.faa") into muscleRefinedMsas

    """
    muscle \
      -in "${fasta}" \
      -out "${fasta.baseName}.faa" \
      -seqtype protein \
      -refine
    """
}

/*
// Construct multiple sequence alignments from the clusters.
process estimateTrees {
    container "quay.io/biocontainers/fasttree:2.1.10--h470a237_2"
    publishDir "trees"

    input:
    set val(type), file(msa) from clustersMSAFastaIndiv.transpose()

    output:
    file "${msa.baseName}.nwk" into indivTrees

    """
    fasttree -fastest -quiet < "${msa}" > "${msa.baseName}.nwk"
    """
}
*/


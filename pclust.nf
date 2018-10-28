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


/*
 * Tidy GFF3s so that genometools doesn't panic.
 */
process tidyGFFs {
    container "quay.io/biocontainers/genometools-genometools:1.5.10--h470a237_1"

    input:
    set val(label), file(fasta), file(gff) from genomes

    output:
    set val(label), file(fasta), file("${gff.baseName}_tidied.gff3") into genomesTidied

    """
    gt gff3 \
      -tidy \
      -sort \
      "${gff}" \
    > "${gff.baseName}_tidied.gff3"
    """
}

genomesTidied.tap { genomes4proteindb; genomes4genomedb }


/*
 * Extract protein sequences from genomes and GFF.
 */
process extractProteins {
    container "quay.io/biocontainers/genometools-genometools:1.5.10--h470a237_1"

    input:
    set val(label), file(fasta), file(gff) from genomes4proteindb

    output:
    file "${label}.faa" optional true into proteins

    """
    gt extractfeat \
      -type CDS \
      -translate \
      -matchdescstart \
      -join \
      -retainids \
      -seqfile "${fasta}" \
      "${gff}" \
    | sed "s/>\\s*/>${label}./g" \
    > "${label}.faa"
    
    if [ ! -s "${label}.faa" ]; then
        rm -f "${label}.faa"
    fi
    """
}


/*
 * Combine all proteins into a single fasta file.
 */
process combineFasta {
    publishDir "sequences"

    input:
    file "*.faa" from proteins.collect()

    output:
    file "combined.faa" into combinedFasta

    """
    cat *faa > combined.faa
    """
}


/*
 * Create the mmseqs2 sequence database
 */
process createSequenceDB {
    container "soedinglab/mmseqs2"
    label 'mmseqs'
    publishDir "sequences"

    input:
    file fasta from combinedFasta

    output:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
        "sequence_h", "sequence_h.index" into seq

    """
    mmseqs createdb "${fasta}" sequence --max-seq-len 14000
    """
}

seq.into {
    seq4Dedup
}


/*
 * Select only unique sequences to cluster.
 */
process getDedupSequences {
    container "soedinglab/mmseqs2"
    label 'mmseqs'
    publishDir "dedup"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
        "sequence_h", "sequence_h.index" from seq4Dedup

    output:
    file "dedup.fasta" into dedupSeqFasta
    file "dedup.tsv" into dedupCluTSV

    """
    mmseqs clusthash sequence result --min-seq-id 1.0
    mmseqs clust sequence result dedup_clu

    mmseqs result2repseq sequence dedup_clu dedup_rep
    mmseqs result2flat sequence sequence dedup_rep dedup.fasta --use-fasta-header

    mmseqs createtsv sequence sequence dedup_clu dedup.tsv
    """
}


/*
 * Create the mmseqs2 sequence database
 */
process createDedupDB {
    container "soedinglab/mmseqs2"
    label 'mmseqs'
    publishDir "dedup"

    input:
    file fasta from dedupSeqFasta

    output:
    set "dedup", "dedup.dbtype", "dedup.index", "dedup.lookup",
        "dedup_h", "dedup_h.index" into dedupSeq

    """
    mmseqs createdb "${fasta}" dedup --max-seq-len 14000
    """
}

dedupSeq.into {
    dedupSeq4Cascade;
    dedupSeq4Profile;
    dedupSeq4Stats;
    dedupSeq4MmseqsMSA;
    dedupSeq4GenomeSearch;
}


/*
 * Perform the first pass clustering using basic mmseqs workflow.
 */
process clusterCascade {
    container "soedinglab/mmseqs2"
    label 'mmseqs'
    publishDir "clusters"

    input:
    set "dedup", "dedup.dbtype", "dedup.index", "dedup.lookup",
        "dedup_h", "dedup_h.index" from dedupSeq4Cascade

    output:
    set "cascade_clusters", "cascade_clusters.index" into cascadeClu

    """
    mkdir -p tmp
    mmseqs cluster \
      dedup \
      cascade_clusters \
      tmp \
      -s 5 \
      --cluster-steps 5 \
      --min-seq-id 0.3
    """
}

cascadeClu.into {
    cascadeClu4Profile;
    cascadeClu4Combine;
}


/*
 * Perform the second clustering pass using sequence profiles.
 * This gets clusters down to about 10%-20% identity.
 */
process clusterProfile {
    container "soedinglab/mmseqs2"
    label 'mmseqs'
    publishDir "clusters"

    input:
    set "cascade_clusters", "cascade_clusters.index" from cascadeClu4Profile
    set "dedup", "dedup.dbtype", "dedup.index", "dedup.lookup",
        "dedup_h", "dedup_h.index" from dedupSeq4Profile

    output:
    set "profile_clusters", "profile_clusters.index" into profileClu

    """
    # Create profiles for each cluster.
    # Generates the profile_consensus file too.
    mmseqs result2profile dedup dedup cascade_clusters profile

    # Search the profiles against the profile consensus sequences.
    # Uses an iterative strategy.
    mkdir -p tmp
    mmseqs search \
      profile \
      profile_consensus \
      result \
      tmp \
      -s 7.5 \
      -c 0.60 \
      -e 0.05 \
      --add-self-matches \
      --num-iterations 3

    # Cluster the matches of profiles vs consensus sequences.
    mmseqs clust profile result profile_clusters_consensus

    # Merge the original clusters with the new ones to get all sequences back.
    mmseqs mergeclusters \
      dedup \
      profile_clusters \
      cascade_clusters \
      profile_clusters_consensus
    """
}

profileClu.into {
    profileClu4Combine;
    profileClu4MSA;
    profileClu4GenomeSearch;
}


/*
 * Merge the two clustering results into a single channel that we can look at.
 */
combinedClu = cascadeClu4Combine.concat(profileClu4Combine)


/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractClusterStats {
    container "soedinglab/mmseqs2"
    label 'mmseqs'
    publishDir { "clusters" }

    input:
    set file(clusters), file(clusters_index) from combinedClu
    set "dedup", "dedup.dbtype", "dedup.index", "dedup.lookup",
        "dedup_h", "dedup_h.index" from dedupSeq4Stats

    output:
    file "${clusters}.tsv" into combinedCluTSV
    file "${clusters}_rep.fasta" into combinedCluRepFasta
    file "${clusters}_stats.tsv" into combinedCluStats

    """
    mmseqs createtsv \
      dedup \
      dedup \
      "${clusters}" \
      "${clusters}.tsv"

    mmseqs result2repseq \
      dedup \
      "${clusters}" \
      "${clusters}_rep"

    mmseqs result2flat \
      dedup \
      dedup \
      "${clusters}_rep" \
      "${clusters}_rep.fasta" \
      --use-fasta-header

    # Do an all vs centroid alignment for each cluster.
    mmseqs align dedup dedup "${clusters}" align -a
    mmseqs convertalis dedup dedup align "${clusters}_stats.tsv"
    """
}


/*
 * Extract sequences from clusters.
 */
process getMmseqsMSA {
    container "soedinglab/mmseqs2"
    label 'mmseqs'

    input:
    set file(clusters), file(clusters_index) from profileClu4MSA
    set "dedup", "dedup.dbtype", "dedup.index", "dedup.lookup",
        "dedup_h", "dedup_h.index" from dedupSeq4MmseqsMSA

    output:
    file "msa.fastalike" into msaFastaLike

    """
    mmseqs result2msa dedup dedup "${clusters}" msa
    mmseqs result2flat dedup dedup msa msa.fastalike
    """
}


/*
 * Extract the individual fasta sequences from the fasta-like file.
 */
process getMmseqsMSAFastas {
    publishDir "msas/mmseqs"

    input:
    file fastalike from msaFastaLike

    output:
    file "*.fasta" into mmseqsMsas

    """
    extract_fastalike.py "${fastalike}"
    """
}


/*
 * Refine the fast MSAs from mmseqs using muscle
 * The issue with the regular mmseqs MSAs is that it can't have gaps in the
 * seed sequence, muscle should refit that.
 */
process refineMSAs {
    container "quay.io/biocontainers/muscle:3.8.1551--h2d50403_3"
    publishDir { "msas/muscle_refine"}

    input:
    file fasta from mmseqsMsas.flatten()

    output:
    file "${fasta.baseName}.faa" into refinedMsas

    """
    NSEQS=\$(grep -c ">" ${fasta})

    if [ "\${NSEQS}" -lt "2" ]; then
      cp "${fasta}" "${fasta.baseName}.faa"
    else
      muscle \
        -in "${fasta}" \
        -out "${fasta.baseName}.faa" \
        -seqtype protein \
        -refine
    fi
    """
}


/*
 * Combine all genomes into a single fasta file.
 */
process combineGenomeFasta {

    input:
    file "*.fasta" from genomes4genomedb.map {l, f, g -> f}.collect()

    output:
    file "combined.fasta" into combinedGenomeFasta

    """
    cat *fasta > combined.fasta
    """
}


/*
 * Construct genomes database
 */
process createGenomeSequenceDB {
    container "soedinglab/mmseqs2"
    label 'mmseqs'

    input:
    file fasta from combinedGenomeFasta

    output:
    set "genome", "genome.dbtype", "genome.index", "genome.lookup",
        "genome_h", "genome_h.index" into genomeSeq

    """
    mmseqs createdb "${fasta}" genome --dont-split-seq-by-len
    """
}


/*
 * Search for the clusters in the original genome sequences.
 */
process searchGenomes {
    container "soedinglab/mmseqs2"
    label 'mmseqs'
    publishDir "genome_search"

    input:
    set file(clusters), file(clusters_index) from profileClu4GenomeSearch
    set "genome", "genome.dbtype", "genome.index", "genome.lookup",
        "genome_h", "genome_h.index" from genomeSeq
    set "dedup", "dedup.dbtype", "dedup.index", "dedup.lookup",
        "dedup_h", "dedup_h.index" from dedupSeq4GenomeSearch

    output:
    file "result.tsv" into searchResults

    """
    mkdir -p tmp
    # Create profiles for each cluster.
    mmseqs result2profile dedup dedup "${clusters}" profile
    # Search profile vs genome
    mmseqs search profile genome result tmp
    # Extract matches from results
    mmseqs convertalis profile genome result result.tsv
    """

}


/*
 * Estimate trees using the MSAs
 */
process estimateTrees {
    container "quay.io/biocontainers/fasttree:2.1.10--h470a237_2"
    publishDir "trees"

    input:
    file msa from refinedMsas

    output:
    file "${msa.baseName}.nwk" into indivTrees

    """
    fasttree -fastest -quiet < "${msa}" > "${msa.baseName}.nwk"
    """
}

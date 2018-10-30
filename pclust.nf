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
    label "genometools"

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
    label "genometools"

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
    label "posix"
    publishDir "sequences"

    input:
    file "*.faa" from proteins.collect()

    output:
    file "proteins.faa" into combinedFasta

    """
    cat *faa > proteins.faa
    """
}


/*
 * Create the mmseqs2 sequence database
 */
process createSequenceDB {
    label 'mmseqs'
    publishDir "sequences"

    input:
    file fasta from combinedFasta

    output:
    file "proteins" into seq

    """
    mkdir -p proteins
    mmseqs createdb "${fasta}" proteins/db --max-seq-len 14000
    """
}

seq.into { seq4Dedup; seq4DedupExtract }


/*
 * Select only unique sequences to cluster.
 */
process clusterDedup {
    label 'mmseqs'
    publishDir "clusters"

    input:
    file "sequence" from seq4Dedup

    output:
    file "dedup" into dedupClu
    file "dedup.tsv" into dedupCluTSV

    """
    mmseqs clusthash sequence/db result --min-seq-id 1.0

    mkdir -p dedup
    mmseqs clust sequence/db result dedup/db

    mmseqs createtsv sequence/db sequence/db dedup/db dedup.tsv
    sed -i '1i cluster\tmember' dedup.tsv
    """
}

dedupClu.set { dedupClu4Extract }


/*
 * Select only unique sequences to cluster.
 */
process getDedupSequences {
    label 'mmseqs'
    publishDir "sequences"

    input:
    file "sequence" from seq4DedupExtract
    file "cluster" from dedupClu4Extract

    output:
    file "dedup.fasta" into dedupSeqFasta

    """
    mmseqs result2repseq sequence/db cluster/db dedup
    mmseqs result2flat sequence/db sequence/db dedup dedup.fasta --use-fasta-header
    """
}


/*
 * Create the mmseqs2 sequence database
 */
process createDedupDB {
    label 'mmseqs'
    publishDir "sequences"

    input:
    file fasta from dedupSeqFasta

    output:
    file "dedup" into dedupSeq

    """
    mkdir -p dedup
    mmseqs createdb "${fasta}" dedup/db --max-seq-len 14000
    """
}

dedupSeq.into {
    dedupSeq4Cascade;
    dedupSeq4Profile;
    dedupSeq4CascadeStats;
    dedupSeq4ProfileStats;
    dedupSeq4MmseqsMSA;
    dedupSeq4GenomeSearch;
}


/*
 * Perform the first pass clustering using basic mmseqs workflow.
 */
process clusterCascade {
    label 'mmseqs'
    publishDir "clusters"

    input:
    file "dedup" from dedupSeq4Cascade

    output:
    file "cascade" into cascadeClu

    """
    mkdir -p cascade

    mkdir -p tmp
    mmseqs cluster \
      dedup/db \
      cascade/db \
      tmp \
      -c 0.8 \
      -s 5 \
      --cluster-steps 4 \
      --min-seq-id 0.3
    """
}

cascadeClu.into {
    cascadeClu4Profile;
    cascadeClu4Stats;
}


/*
 * Perform the second clustering pass using sequence profiles.
 * This gets clusters down to about 10%-20% identity.
 */
process clusterProfile {
    label 'mmseqs'
    publishDir "clusters"

    input:
    file "cascade" from cascadeClu4Profile
    file "dedup" from dedupSeq4Profile

    output:
    file "profile" into profileClu

    """
    # Create profiles for each cluster.
    # Generates the profile_consensus file too.
    mmseqs result2profile dedup/db dedup/db cascade/db profile1

    # Search the profiles against the profile consensus sequences.
    # Uses an iterative strategy.
    mkdir -p tmp
    mmseqs search \
      profile1 \
      profile1_consensus \
      result \
      tmp \
      -s 7.5 \
      -c 0.80 \
      -e 0.05 \
      --add-self-matches \
      --num-iterations 3

    # Cluster the matches of profiles vs consensus sequences.
    mmseqs clust profile1 result profile1_clusters_consensus

    # Merge the original clusters with the new ones to get all sequences back.
    mkdir -p profile
    mmseqs mergeclusters \
      dedup/db \
      profile/db \
      cascade/db \
      profile1_clusters_consensus
    """
}

profileClu.into {
    profileClu4Stats;
    profileClu4MSA;
    profileClu4GenomeSearch;
}


/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractCascadeClusterStats {
    label 'mmseqs'
    publishDir "clusters"

    input:
    file clusters from cascadeClu4Stats
    file seq from dedupSeq4CascadeStats

    output:
    file "${clusters}.tsv" into cascadeCluTSV
    file "${clusters}_rep.fasta" into cascadeCluRepFasta
    file "${clusters}_stats.tsv" into cascadeCluStats

    script:
    template "mmseqs_cluster_stats.sh"
}


/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractClusterStats {
    label 'mmseqs'
    publishDir "clusters"

    input:
    file clusters from profileClu4Stats
    file seq from dedupSeq4ProfileStats

    output:
    file "${clusters}.tsv" into profileCluTSV
    file "${clusters}_rep.fasta" into profileCluRepFasta
    file "${clusters}_stats.tsv" into profileCluStats

    script:
    template "mmseqs_cluster_stats.sh"
}


/*
 * Extract sequences from clusters.
 */
process getMmseqsMSA {
    label 'mmseqs'

    input:
    file "profile" from profileClu4MSA
    file "dedup" from dedupSeq4MmseqsMSA

    output:
    file "msa.fastalike" into msaFastaLike

    """
    mmseqs result2msa dedup/db dedup/db profile/db msa
    mmseqs result2flat dedup/db dedup/db msa msa.fastalike
    """
}


/*
 * Extract the individual fasta sequences from the fasta-like file.
 */
process getMmseqsMSAFastas {
    label "python3"
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
    label "muscle"
    publishDir "msas/muscle"

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
 * Estimate trees using the MSAs
 */
process estimateTrees {
    label "fasttree"
    publishDir "msas/trees"

    input:
    file msa from refinedMsas

    output:
    file "${msa.baseName}.nwk" into indivTrees

    """
    fasttree -fastest -quiet < "${msa}" > "${msa.baseName}.nwk"
    """
}


/*
 * Add the genome name to the scaffold names in the genome fasta.
 */
process addGenomeNameToScaffold {
    label "posix"
 
    input:
    set val(label), file(fasta), file(gff) from genomes4genomedb

    output:
    file "genome.fasta" into genomesWithNames

    """
    LABEL=\$(basename ${label})
    sed "s/>\\s*/>\${LABEL}./g" < ${fasta} > genome.fasta
    """
}


/*
 * Combine all genomes into a single fasta file.
 */
process combineGenomeFasta {
    label "posix"
    publishDir "sequences"

    input:
    file "*.fasta" from genomesWithNames.collect()

    output:
    file "genomes.fasta" into combinedGenomeFasta

    """
    cat *.fasta > genomes.fasta
    """
}


/*
 * Construct genomes database
 */
process createGenomeSequenceDB {
    label 'mmseqs'
    publishDir "sequences"

    input:
    file fasta from combinedGenomeFasta

    output:
    file "genomes" into genomeSeq

    """
    mkdir -p genomes
    mmseqs createdb "${fasta}" genomes/db --dont-split-seq-by-len
    """
}


/*
 * Search for the clusters in the original genome sequences.
 */
process searchGenomes {
    label 'mmseqs'
    publishDir "genome_searches"

    input:
    file "cascade" from profileClu4GenomeSearch
    file "genomes" from genomeSeq
    file "dedup" from dedupSeq4GenomeSearch

    output:
    file "profile_matches.tsv" into searchResults

    """
    mkdir -p tmp
    # Create profiles for each cluster.
    mmseqs result2profile dedup/db dedup/db cascade/db profile

    # Search profile vs genome
    # Search parameters are slightly more conservative than default.
    mmseqs search \
      profile \
      genomes/db \
      result \
      tmp \
      --realign \
      --gapopen 15 \
      --gapextend 2 \
      --cov-mode 2 \
      --rescore-mode 2 \
      --min-length 20 \
      --orf-start-mode 1 \
      --use-all-table-starts true

    # Extract matches from results
    mmseqs convertalis \
      profile \
      genomes/db \
      result profile_matches.tsv \
      --format-mode 0 \
      --format-output "query target evalue qcov tcov gapopen pident nident mismatch raw bits qstart qend tstart tend qlen tlen alnlen cigar qframe tframe"

    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen\tcigar\tqframe\ttframe' profile_matches.tsv
    """
}

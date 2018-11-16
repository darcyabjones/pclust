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

params.trees = false
params.nomsa = false
params.global_profile = false
params.global_clusters = false
params.global_seqs = false

proteins = Channel.fromPath( params.proteins )

if (params.global_profile) {
    globalProfile = Channel.fromPath( params.global_profile )
} else if ( params.global_clusters && params.global_seqs ) {
    globalClusters = Channel.fromPath( params.global_clusters )
    globalSeqs = Channel.fromPath( params.global_seqs )
} else {
    log.info "Need either the global profile or global clusters + seqdb"
    exit 1
}


if ( params.trees && params.nomsa ) {
    log.info "Cannot construct trees without msas"
    exit 1
}


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


seq.into {
    seq4Cascade;
    seq4CreateLocalProfile;
    seq4Profile;
    seq4CascadeStats;
    seq4ProfileStats;
    seq4MmseqsMSA;
    seq4SearchGlobal
}


/*
 * Perform the first pass clustering using basic mmseqs workflow.
 */
process clusterCascade {
    label 'mmseqs'
    publishDir "clusters"

    input:
    file "seq" from seq4Cascade

    output:
    file "cascade" into cascadeClu

    """
    mkdir -p cascade

    mkdir -p tmp
    mmseqs cluster \
      seq/db \
      cascade/db \
      tmp \
      --threads ${task.cpus} \
      -c 0.7 \
      --cov-mode 0 \
      -s 6 \
      --cluster-steps 3 \
      --cluster-mode 0 \
      --min-seq-id 0.3
    """
}

cascadeClu.into {
    cascadeClu4CreateLocalProfile;
    cascadeClu4Profile;
    cascadeClu4Stats;
    cascadeClu4MSA;
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
    file seq from seq4CascadeStats

    output:
    file "${clusters}.tsv" into cascadeCluTSV
    file "${clusters}_rep.fasta" into cascadeCluRepFasta
    file "${clusters}_stats.tsv" into cascadeCluStats

    script:
    template "mmseqs_cluster_stats.sh"
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

globalProfile.into { globalProfile4Search; globalProfile4ExtractProfile }

process searchGlobal {
    label "mmseqs"
    publishDir "clusters"

    input:
    file "local_seqs" from seq4SearchGlobal
    file "global_profile" from globalProfile4Search

    output:
    file "search_global" into globalSearchResults 
    file "search_global.tsv" into globalSearchResultsTsv

    """
    mkdir -p tmp
    mkdir -p search_global
    mmseqs search local_seqs/db global_profile/db search_global/db tmp

    mmseqs convertalis \
      local_seqs/db \
      global_profile/db \
      search_global/db \
      "search_global.tsv" \
      --threads ${task.cpus} \
      --format-mode 0 \
      --format-output "query target evalue qcov tcov gapopen pident nident mismatch raw bits qstart qend tstart tend qlen tlen alnlen cigar"

    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen\tcigar' search_global.tsv
    """
}


process createLocalProfile {
    label "mmseqs"
    publishDir "clusters"

    input:
    file "local_seqs" from seq4CreateLocalProfile
    file "global_profile" from globalProfile4ExtractProfile
    file "search_global" from globalSearchResults

    output:
    file "local_profile" into localProfile   
 
    """
    mkdir -p local_profile
    mmseqs result2profile local_seqs/db global_profile/db search_global/db local_profile/db
    """
}


/*
 * Perform the second clustering pass using sequence profiles.
 * This gets clusters down to about 10%-20% identity.
 */
process clusterProfile {
    label 'mmseqs'
    publishDir "clusters"

    input:
    file "local_profile" from localProfile

    output:
    file "profile" into profileClu
    file "profile_matches.tsv" into profileSearchResults

    """
    # Search the profiles against the profile consensus sequences.
    # Uses an iterative strategy.
    mkdir -p tmp
    mkdir search_result
    mmseqs search \
      local_profile/db \
      local_profile/db_consensus \
      search_result/db \
      tmp \
      --threads ${task.cpus} \
      -s 7.0 \
      -c 0.60 \
      --add-self-matches \
      --num-iterations 3

    # Cluster the matches of profiles vs consensus sequences.
    mkdir -p profile
    mmseqs clust \
      local_profile/db \
      search_result/db \
      profile/db \
      --threads ${task.cpus}

    mmseqs convertalis \
      local_profile/db \
      local_profile/db \
      search_result/db \
      "profile_matches.tsv" \
      --threads ${task.cpus} \
      --format-mode 0 \
      --format-output "query target evalue qcov tcov gapopen pident nident mismatch raw bits qstart qend tstart tend qlen tlen alnlen cigar"

    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen\tcigar' profile_matches.tsv
    """
}


profileClu.into {
    profileClu4Stats;
    profileClu4MSA;
}

/
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractProfileClusterStats {
    label 'mmseqs'
    publishDir "clusters"

    input:
    file clusters from profileClu4Stats
    file seq from seq4ProfileStats

    output:
    file "${clusters}.tsv" into profileCluTSV
    file "${clusters}_rep.fasta" into profileCluRepFasta
    file "${clusters}_stats.tsv" into profileCluStats

    script:
    template "mmseqs_cluster_stats.sh"
}


/*
 * Merge the cluster tables and profile statistics
 */
process joinClusterStats {
    label 'R'
    publishDir "clusters"

    input:
    file "cascade.tsv" from cascadeCluTSV
    file "profile.tsv" from profileCluTSV
    file "profile_stats.tsv" from profileCluStats

    output:
    file "clusters.tsv" into clusterStats

    """
    join_clusters.R \
      cascade.tsv \
      profile.tsv \
      profile_stats.tsv \
    > clusters.tsv
    """
}


if ( !params.nomsa ) {
    /*
     * Extract sequences from clusters.
     */
    process getMmseqsMSA {
        label 'mmseqs'

        input:
        file "clusters" from clu4MSAs
        file "seq" from seq4MmseqsMSA

        output:
        file "msa.fastalike" into msaFastaLike

        """
        mmseqs result2msa \
          seq/db \
          seq/db \
          clusters/db \
          msa \
          --threads ${task.cpus}

        mmseqs result2flat \
          seq/db \
          seq/db \
          msa \
          msa.fastalike
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
        file "*.fasta" into mmseqsMsas4Refinement

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
        tag { fasta.baseName }

        input:
        file fasta from mmseqsMsas4Refinement

        output:
        file "${fasta.baseName}.faa" into refinedMsas

        """
        muscle \
          -in "${fasta}" \
          -out "${fasta.baseName}.faa" \
          -seqtype protein \
          -refine
        """
    }
    
    msas4Trees = refinedMsas
}


if ( params.trees ) {
    /*
     * Estimate trees using the MSAs
     */
    process estimateTrees {
        label "fasttree"
        publishDir "msas/trees"

        input:
        file msa from msas4Trees

        output:
        file "${msa.baseName}.nwk" into indivTrees

        """
        FastTree -fastest -quiet < "${msa}" > "${msa.baseName}.nwk"
        """
    }
}

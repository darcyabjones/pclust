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
      --proteins               description

    Options:
      --profile
      --nomsa          description
      --nomsa_refine          description
      --tree         description

    Outputs:

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}


params.profile = false
params.nomsa = false
params.nomsa_refine = false
params.tree = false


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


seq.into {
    seq4Cascade;
    seq4Profile;
    seq4CascadeStats;
    seq4ProfileStats;
    seq4MmseqsMSA;
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
      -c 0.8 \
      --cov-mode 0 \
      -s 6 \
      --cluster-steps 4 \
      --cluster-mode 0 \
      --min-seq-id 0.1
    """
}

cascadeClu.into {
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


if ( params.profile ) {
    /*
     * Perform the second clustering pass using sequence profiles.
     * This gets clusters down to about 10%-20% identity.
     */
    process clusterProfile {
        label 'mmseqs'
        publishDir "clusters"

        input:
        file "cascade" from cascadeClu4Profile
        file "seq" from seq4Profile

        output:
        file "profile" into profileClu

        """
        # Create profiles for each cluster.
        # Generates the profile_consensus file too.
        mmseqs result2profile \
          seq/db \
          seq/db \
          cascade/db \
          profile1 \
          --threads ${task.cpus}

        # Search the profiles against the profile consensus sequences.
        # Uses an iterative strategy.
        mkdir -p tmp
        mmseqs search \
          profile1 \
          profile1_consensus \
          result \
          tmp \
          --threads ${task.cpus} \
          -s 7.0 \
          -c 0.70 \
          --add-self-matches \
          --num-iterations 2

        # Cluster the matches of profiles vs consensus sequences.
        mmseqs clust \
          profile1 \
          result \
          profile1_clusters_consensus \
          --threads ${task.cpus}

        # Merge the original clusters with the new ones to get all sequences back.
        mkdir -p profile
        mmseqs mergeclusters \
          seq/db \
          profile/db \
          cascade/db \
          profile1_clusters_consensus \
          --threads ${task.cpus}
        """
    }

    profileClu.into {
        profileClu4Stats;
        profileClu4MSA;
    }


    /*
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

    clu4MSAs = profileClu4MSA
} else {
    clu4MSAs = cascadeClu4MSA
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
        file "*.fasta" into mmseqsMsas

        """
        extract_fastalike.py "${fastalike}"
        """
    }
}

if ( !params.nomsa_refine ) {
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
}

if ( params.tree && !params.nomsa_refine ) {
    msas4Trees = refinedMsas
} else if ( params.tree && !params.nomsa ) {
    msas4Trees = mmseqsMsas
}

if ( params.tree ) {
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

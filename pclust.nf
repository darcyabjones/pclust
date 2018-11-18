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
      --seqs              description
      --db
      --enrich_seqs
      --enrich_db

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
params.nomsa_refine = false
params.enrich_db = false
params.enrich_seqs = false
params.seqs = false
params.db = false


if ( params.trees && params.nomsa ) {
    log.info "Cannot construct trees without msas"
    exit 1
}


if ( params.db ) {
    seqDB = Channel.fromPath( params.db )
} else if ( params.seqs ) {
    proteins = Channel.fromPath( params.seqs )

    /*
     * Create the mmseqs2 sequence database
     */
    process createSequenceDB {
        label 'mmseqs'
        publishDir "sequences"
    
        input:
        file fasta from proteins
    
        output:
        file "seqdb" into seqDB
    
        """
        mkdir -p seqdb
        mmseqs createdb "${fasta}" seqdb/db --max-seq-len 14000
        """
    }
} else {
    log.info "Please provide either sequences or the seqdb"
    exit 1
}


seqDB.into {
    seqDB4ClusterHighId;
    seqDB4CreateHighIdSubDB;
    seqDB4MergeClusters;
    seqDB4CascadeStats;
    seqDB4CreateProfile;
    seqDB4ClusterProfile;
    seqDB4ProfileStats;
    seqDB4MmseqsMSA;
}


if (params.enrich_db) {
    enrichSeqsDB = Channel.fromPath( params.enrich_db )
    enrich = true
} else if ( params.enrich_seqs ) {
    enrichSeqs = Channel.fromPath( params.enrich_seqs )
    enrich = true

    process createEnrichSeqsDB {
        label 'mmseqs'
        publishDir "sequences"

        input:
        file fasta from enrichSeqs

        output:
        file "enrich_db" into enrichSeqsDB

        """
        mkdir -p enrich_db
        mmseqs createdb "${fasta}" enrich_db/db --max-seq-len 14000
        """
    }
} else {
    enrich = false
}


process clusterHighId {
    label 'mmseqs'
    publishDir "clusters"

    input:
    file "seq" from seqDB4ClusterHighId

    output:
    file "high_id" into highIdClu

    """
    mkdir -p "high_id"
    mkdir -p "tmp"

    mmseqs linclust \
      seq/db \
      high_id/db \
      tmp \
      --min-seq-id 0.90 \
      -c 0.8 \
      --cov-mode 0

    rm -rf -- tmp
    """
}

highIdClu.into {
    highIdClu4CreateHighIdSubDB;
    highIdClu4MergeClusters;
}


process createHighIdSubDB {
    label 'mmseqs'

    input:
    file "high_id" from highIdClu4CreateHighIdSubDB
    file "seq" from seqDB4CreateHighIdSubDB

    output:
    file "high_id_db" into highIdSubDB

    """
    mkdir -p high_id_db
    mmseqs createsubdb "high_id/db" "seq/db" "high_id_db/db"
    """
}

highIdSubDB.into {
    highIdSubDB4ClusterCascade;
}


/*
 * Perform the first pass clustering using basic mmseqs workflow.
 */
process clusterCascade {
    label 'mmseqs'

    input:
    file "seq" from highIdSubDB4ClusterCascade

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
      --min-seq-id 0.3 \
      -c 0.7 \
      --cov-mode 0 \
      -s 7.5 \
      --cluster-steps 5 \
      --cluster-mode 0 \
      --max-seqs 100

    rm -rf -- tmp
    """
}

cascadeClu.set { cascadeClu4MergeClusters }

process mergeClusters {
    label "mmseqs"
    publishDir "clusters"

    input:
    file "seq" from seqDB4MergeClusters
    file "linclust" from highIdClu4MergeClusters
    file "cascade_clu" from cascadeClu4MergeClusters

    output:
    file "cascade" into mergedClu

    """
    mkdir -p cascade

    mmseqs mergeclusters \
      seq/db \
      cascade/db \
      linclust/db \
      cascade_clu/db
    """
}


mergedClu.into {
    cascadeClu4CreateProfile;
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
    file seq from seqDB4CascadeStats

    output:
    file "${clusters}.tsv" into cascadeCluTSV
    file "${clusters}_rep.fasta" into cascadeCluRepFasta
    file "${clusters}_stats.tsv" into cascadeCluStats

    script:
    template "mmseqs_cluster_stats.sh"
}
    

process createProfile {
    label "mmseqs"
    publishDir "clusters"
    
    input:
    file "clusters" from cascadeClu4CreateProfile
    file "seqs" from seqDB4CreateProfile
    
    output:
    file "profile" into profile
    
    """
    mkdir -p profile
    
    # Create profiles for each cluster.
    # Generates the profile_consensus file too.
    mmseqs result2profile \
      seqs/db \
      seqs/db \
      clusters/db \
      profile/db \
      --threads ${task.cpus}
    """
}


if (! enrich ) {
    profile.set { profile4Clu }
} else {
    profile.into { profile4Search; profile4EnrichProfile }
    enrichSeqsDB.into { enrichSeqsDB4Search; enrichSeqsDB4EnrichProfile }

    
    process enrichProfile {
        label "mmseqs"
    
        input:
        file "profile" from profile4Search
        file "enrich_seqs" from enrichSeqsDB4Search
    
        output:
        file "enrich_matches" into enrichSearchResults 
    
        """
        mkdir -p tmp
        mkdir -p enrich_matches
        mmseqs search \
          profile/db \
          enrich_seqs/db \
          enrich_matches/db \
          tmp \
          --max-seqs 1000 \
          -a \
          -e 0.00001 \
          --e-profile 0.01 \
          -c 0.1 \
          -s 7.5 \
          --rescore-mode 1 \
          --num-iterations 3
        """
    }
    

    process createEnrichedProfile {
        label "mmseqs"
    
        input:
        file "input_profile" from profile4EnrichProfile
        file "enrich_seqs" from enrichSeqsDB4EnrichProfile
        file "enrich_matches" from enrichSearchResults
    
        output:
        file "enriched_profile" into enrichedProfile   
     
        """
        mkdir -p enriched_profile
        mmseqs result2profile \
          input_profile/db \
          enrich_seqs/db \
          enrich_matches/db \
          enriched_profile/db
        """
    }
    enrichedProfile.set { profile4Clu }
} 


/*
 * Perform the second clustering pass using sequence profiles.
 * This gets clusters down to about 10%-20% identity.
 */
process clusterProfile {
    label 'mmseqs'
    publishDir "clusters"

    input:
    file "input_profile" from profile4Clu

    output:
    file "profile" into profileClu
    file "profile_matches.tsv" into profileSearchResults

    """
    # Search the profiles against the profile consensus sequences.
    # Uses an iterative strategy.
    mkdir -p tmp
    mkdir search_result
    mmseqs search \
      input_profile/db \
      input_profile/db_consensus \
      search_result/db \
      tmp \
      --threads ${task.cpus} \
      --max-seqs 1000 \
      -c 0.6 \
      --cov-mode 0 \
      -s 7.5 \
      --cov 0.6 \
      --add-self-matches \
      --num-iterations 2
      

    # Cluster the matches of profiles vs consensus sequences.
    mkdir -p profile
    mmseqs clust \
      input_profile/db \
      search_result/db \
      profile/db \
      --threads ${task.cpus} \
      --cluster-mode 1

    mmseqs convertalis \
      input_profile/db \
      input_profile/db_consensus \
      search_result/db \
      "profile_matches.tsv" \
      --threads ${task.cpus} \
      --format-mode 0 \
      --format-output "query target evalue qcov tcov gapopen pident nident mismatch raw bits qstart qend tstart tend qlen tlen alnlen"

    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' profile_matches.tsv
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
    file seq from seqDB4ProfileStats

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
        file "clusters" from profileClu4MSA
        file "seq" from seqDB4MmseqsMSA

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
}

if ( !params.nomsa_refine && !params.nomsa ) {
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
        file fasta from mmseqsMsas4Refinement.flatten()

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

    msas4Trees = refinedMsas
} else if ( !params.nomsa ) {
    msas4Trees = mmseqsMsas4Refinement.flatten()
}


if ( params.trees ) {
    /*
     * Estimate trees using the MSAs
     */
    process estimateTrees {
        label "fasttree"
        publishDir "msas/trees"
        tag { msa.baseName }

        input:
        file msa from msas4Trees

        output:
        file "${msa.baseName}.nwk" into indivTrees

        """
        FastTree -fastest -quiet < "${msa}" > "${msa.baseName}.nwk"
        """
    }
}

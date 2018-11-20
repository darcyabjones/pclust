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
      --seqs              description
      --db
      --enrich_seqs
      --enrich_db

    Options:
      --trees
      --nomsa
      --nomsa_refine
      --enrich_seqs
      --enrich_db
      --enrich_msa

    Outputs:

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

params.seqs = false
params.db = false
params.trees = false
params.nomsa = false
params.nomsa_refine = false
params.enrich_db = false
params.enrich_seqs = false
params.enrich_profile = false
params.enrich_msa = false


if ( params.trees && params.nomsa ) {
    log.info "Cannot construct trees without msas"
    exit 1
}

if ( params.enrich_msa && !(params.enrich_seqs || params.enrich_db)) {
    log.info "Enriching msas requires a database or the sequences"
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
        publishDir "${params.outdir}/sequences"
    
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
    seqDB4ClusterCascade;
    seqDB4MergeClusters;
    seqDB4CascadeStats;
    seqDB4CreateProfile;
    seqDB4ClusterProfile;
    seqDB4ProfileStats;
    seqDB4CreateProfileCluProfile;   
    seqDB4MmseqsMSA;
}


if (params.enrich_db && ( params.enrich_profile || params.enrich_msa )) {
    enrichSeqsDB = Channel.fromPath( params.enrich_db )
} else if ( params.enrich_seqs && ( params.enrich_profile || params.enrich_msa )) {
    enrichSeqs = Channel.fromPath( params.enrich_seqs )

    process createEnrichSeqsDB {
        label 'mmseqs'
        publishDir "${params.outdir}/sequences"

        input:
        file fasta from enrichSeqs

        output:
        file "enrich_db" into enrichSeqsDB

        """
        mkdir -p enrich_db
        mmseqs createdb "${fasta}" enrich_db/db --max-seq-len 14000
        """
    }
} else if ( params.enrich_profile || params.enrich_msa ) {
    log.info "You asked to enrich the profile and/or the msa but didn't provide db to enrich with."
    exit 1
}


/*
 * Perform fast high-identity clustering.
 */
process clusterHighId {
    label 'mmseqs'
    publishDir "${params.outdir}/clusters"

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


/*
 * Extract representative sequences of high-identity clustering as seq database.
 */
process createHighIdSubDB {
    label 'mmseqs'
    publishDir "${params.outdir}/clusters"

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

highIdSubDB.set {
    highIdSubDB4ClusterCascade;
}


/*
 * Perform the first pass clustering using basic mmseqs workflow.
 */
process clusterCascade {
    label 'mmseqs'

    input:
    file "seq" from seqDB4ClusterCascade //highIdSubDB4ClusterCascade

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
      --min-seq-id 0.0 \
      -c 0.7 \
      --cov-mode 0 \
      -s 7.5 \
      --cluster-steps 5 \
      --cluster-mode 0

    rm -rf -- tmp
    """
}

cascadeClu.into {
    cascadeClu4CreateProfile;
    cascadeClu4Profile;
    cascadeClu4Stats;
    cascadeClu4MSA;
    cascadeClu4MergeClusters;
}



/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractCascadeClusterStats {
    label 'mmseqs'
    publishDir "${params.outdir}/clusters"

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
    publishDir "${params.outdir}/clusters"
    
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


if ( params.enrich_profile || params.enrich_msa ) {
    enrichSeqsDB.into {
        enrichSeqsDB4Search;
        enrichSeqsDB4EnrichProfile;
        enrichSeqsDB4enrichMSA;
    }
}

if ( params.enrich_profile ) {
    profile.into { profile4Search; profile4EnrichProfile }

    /*
     * Enrich the sequences by searching a database.
     */
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
          --max-seqs 200 \
          -e 0.00001 \
          --e-profile 0.01 \
          -c 0.2 \
          --start-sens 5.0 \
          --sens-steps 2 \
          -s 7.0 \
          --rescore-mode 1 \
          --num-iterations 2

        rm -rf -- tmp
        """
    }
    

    /*
     * Convert search results into an enriched profile.
     */
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
} else {
    profile.set { profile4Clu }
} 


/*
 * Perform the third clustering pass using sequence profiles.
 * This gets clusters down to about 10%-20% identity.
 */
process clusterProfile {
    label 'mmseqs'
    publishDir "${params.outdir}/clusters"

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
      --max-seqs 300 \
      -c 0.8 \
      --cov-mode 1 \
      --start-sens 5 \
      --sens-steps 2 \
      -s 7.0 \
      --add-self-matches \
      --num-iterations 2
      
    # Cluster the matches of profiles vs consensus sequences.
    mkdir -p profile
    mmseqs clust \
      input_profile/db \
      search_result/db \
      profile/db \
      --threads ${task.cpus} \
      --cluster-mode 2

    mmseqs convertalis \
      input_profile/db \
      input_profile/db_consensus \
      search_result/db \
      "profile_matches.tsv" \
      --threads ${task.cpus} \
      --format-mode 0 \
      --format-output "query target evalue qcov tcov gapopen pident nident mismatch raw bits qstart qend tstart tend qlen tlen alnlen"

    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' profile_matches.tsv

    rm -rf -- tmp
    """
}


/*
 * Merge the clustering results into single db.
 */
process mergeClusters {
    label "mmseqs"
    publishDir "${params.outdir}/clusters"

    input:
    file "seq" from seqDB4MergeClusters
    file "cascade" from cascadeClu4MergeClusters
    file "profile" from profileClu

    output:
    file "merged" into mergedClu

    """
    mkdir -p merged

    mmseqs mergeclusters \
      seq/db \
      merged/db \
      cascade/db \
      profile/db
    """
}

mergedClu.into {
    profileClu4Stats;
    profileClu4MSA;
}

/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractProfileClusterStats {
    label 'mmseqs'
    publishDir "${params.outdir}/clusters"

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
process joinClusterStats {
    label 'R'
    publishDir "${params.outdir}/clusters"

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
 */


if ( ( !params.nomsa ) && params.enrich_msa ) {

    /*
     * Create a profile database from the final clusters.
     */
    process createProfileCluProfile {
        label "mmseqs"
        
        input:
        file "clusters" from profileClu4MSA
        file "seqs" from seqDB4CreateProfileCluProfile
        
        output:
        file "profile" into profileCluProfile
        
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

    profileCluProfile.into { profile4EnrichMSAs; profile4EnrichMSAsResults }

    /*
     * Search the enrichment database to enrich msas.
     */
    process enrichMSAs {
        label "mmseqs"
    
        input:
        file "profile" from profile4EnrichMSAs
        file "enrich_seqs" from enrichSeqsDB4enrichMSA
    
        output:
        file "enrich_matches" into enrichMsaSearchResults 
    
        """
        mkdir -p tmp
        mkdir -p enrich_matches
        mmseqs search \
          profile/db \
          enrich_seqs/db \
          enrich_matches/db \
          tmp \
          --max-seqs 1000 \
          -e 0.00001 \
          --e-profile 0.01 \
          -c 0.1 \
          --start-sens 5 \
          --sens-steps 1 \
          -s 7.0 \
          --rescore-mode 1 \
          --num-iterations 3
        rm -rf -- tmp
        """
    }

    enrichMsaSearchResults.into { clu4MSA; enrichMsaSearchResults4Results }
    
    process enrichMSAsResults {
        label "mmseqs"
        publishDir "msas"

        input:
        file "profile" from profile4EnrichMSAsResults
        file "matches" from enrichMsaSearchResults4Results

        output:
        file "enrich_matches.tsv" into enrichedMatches

        """
        mmseqs convertalis \
          input_profile/db \
          input_profile/db_consensus \
          search_result/db \
          "enrich_matches.tsv" \
          --threads ${task.cpus} \
          --format-mode 0 \
          --format-output "query target evalue qcov tcov gapopen pident nident mismatch raw bits qstart qend tstart tend qlen tlen alnlen"

        sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' enrich_matches.tsv
        """
    }

} else {
    clu4MSA = profileClu4MSA
}


if ( !params.nomsa ) {
    /*
     * Extract sequences from clusters/profiles.
     */
    process getMmseqsMSA {
        label 'mmseqs'

        input:
        file "clusters" from clu4MSA
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
        publishDir "${params.outdir}/msas/mmseqs"

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
        publishDir "${params.outdir}/msas/muscle"
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
        publishDir "${params.outdir}/msas/trees"
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

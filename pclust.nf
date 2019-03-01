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

    seqdb = Channel
        .fromPath( params.db, type: 'dir', checkIfExists: true, glob: false )
        .first()

} else if ( params.seqs ) {

    proteins = Channel
        .fromPath( params.seqs, type: 'file', checkIfExists: true, glob: false )
        .first()

    /*
     * Create the mmseqs2 sequence database
     */
    process createSequenceDB {
        label 'mmseqs'
        publishDir "${params.outdir}/sequences"

        input:
        file "seqs.fasta" from proteins

        output:
        file "seqdb" into seqdb

        script:
        FASTA = "seqs.fasta"
        OUTDB = "seqdb"
        template: mmseqs_createdb.sh
    }

} else {
    log.info "Please provide either sequences or the seqdb"
    exit 1
}


if (params.enrich_db && ( params.enrich_profile || params.enrich_msa )) {

    enrichdb = Channel
        .fromPath( params.enrich_db, type: 'dir', checkIfExists: true, glob: false )
        .first()

} else if ( params.enrich_seqs && ( params.enrich_profile || params.enrich_msa )) {

    enrichSeqs = Channel
        .fromPath( params.enrich_seqs, type: 'file', checkIfExists: true, glob: false )
        .first()


    process createEnrichSeqsDB {
        label 'mmseqs'
        publishDir "${params.outdir}/sequences"

        input:
        file "seqs.fasta" from enrichSeqs

        output:
        file "enrich_db" into enrichdb

        script:
        FASTA = "seqs.fasta"
        OUTDB = "enrich_db"
        template 'mmseqs_createdb.sh'
    }

} else if ( params.enrich_profile || params.enrich_msa ) {
    log.info "You asked to enrich the profile and/or the msa but didn't provide db to enrich with."
    exit 1
}


/*
 * Perform the first pass clustering using basic mmseqs workflow.
 */
process clusterCascade {
    label 'mmseqs'

    input:
    file "seq" from seqdb

    output:
    file "cascade" into cascadeClusters

    script:
    INDB = "seq"
    OUTDB = "cascade"
    NCPUS = task.cpus
    template mmseqs_cluster_cascade.sh
}


/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractCascadeClusterStats {
    label 'mmseqs'
    publishDir "${params.outdir}/clusters"

    input:
    file "seqs" from seqdb
    file "clusters" from cascadeClusters

    output:
    file "${clusters}.tsv" into cascadeCluTSV
    file "${clusters}_rep.fasta" into cascadeCluRepFasta
    file "${clusters}_stats.tsv" into cascadeCluStats

    script:
    SEQS = "seqs"
    CLUSTERS = "clusters"
    NCPUS = task.cpus
    template "mmseqs_cluster_stats.sh"
}


process createProfile {
    label "mmseqs"
    publishDir "${params.outdir}/clusters"

    input:
    file "clusters" from cascadeClu4CreateProfile
    file "seqs" from seqdb

    output:
    file "profile" into profile

    script:
    QUERY = "seqs"
    TARGET = "seqs"
    RESULTS = "clusters"
    OUTDB = "profile"
    NCPUS = task.cpus
    template mmseqs_result_to_profiles.sh
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

        script:
        PROFILE = "profile"
        TARGET = "enrich_seqs"
        OUTDB = "enrich_matches"
        NCPUS = task.cpus
        template mmseqs_search_profile_strict.sh
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

        script:
        QUERY = "input_profile"
        TARGET = "enrich_seqs"
        RESULTS = "enrich_matches"
        OUTDB = "enriched_profiles"
        NCPUS = task.cpus
        template mmseqs_result_to_profiles.sh
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

    script:
    INDB = "input_profile"
    OUTDB = "profile"
    NCPUS = task.cpus
    template mmseqs_cluster_profiles.sh
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
    file "clusters" from profileClu4Stats
    file "seq" from seqDB4ProfileStats

    output:
    file "${clusters}.tsv" into profileCluTSV
    file "${clusters}_rep.fasta" into profileCluRepFasta
    file "${clusters}_stats.tsv" into profileCluStats

    script:
    SEQS = "seq"
    CLUSTERS = "clusters"
    NCPUS = task.cpus
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

        script:
        QUERY = "seqs"
        TARGET = "seqs"
        RESULTS = "clusters"
        OUTDB = "profile"
        template mmseqs_result_to_profile.sh
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

        script:
        PROFILE = "profile"
        TARGET = "enrich_seqs"
        OUTDB = "enrich_matches"
        NCPUS = task.cpus
        template: mmseqs_search_profile_relaxed.sh
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



/*

if ( !params.nomsa_refine && !params.nomsa ) {
     * Refine the fast MSAs from mmseqs using muscle
     * The issue with the regular mmseqs MSAs is that it can't have gaps in the
     * seed sequence, muscle should refit that.
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
*/


/*
if ( params.trees ) {
     * Estimate trees using the MSAs
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
*/

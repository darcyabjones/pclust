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
params.enrich_db = false
params.enrich_seqs = false
params.enrich_profile = false
params.enrich_msa = false
params.mafft = false
params.hhdata = false

if ( params.enrich_msa && !(params.enrich_seqs || params.enrich_db)) {
    log.info "Enriching msas requires a database or the sequences"
    exit 1
}

if ( params.enrich_profile && !(params.enrich_seqs || params.enrich_db)) {
    log.info "Enriching profiles for clustering requires a database or the sequences"
    exit 1
}

if ( params.hhdata ) {
    hhdata = Channel
        .fromPath( params.hhdata, type: 'dir', checkIfExists: true, glob: false )
        .first()
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
        template "mmseqs_createdb.sh"
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
        template "mmseqs_createdb.sh"
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
    template "mmseqs_cluster_cascade.sh"
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
    file "cascade" from cascadeClusters

    output:
    file "cascade.tsv" into cascadeClustersTSV
    file "cascade_rep.fasta" into cascadeClustersRepFasta
    file "cascade_stats.tsv" into cascadeClustersStats

    script:
    SEQS = "seqs"
    CLUSTERS = "cascade"
    NCPUS = task.cpus
    template "mmseqs_cluster_stats.sh"
}


process createProfile {
    label "mmseqs"
    publishDir "${params.outdir}/clusters"

    input:
    file "clusters" from cascadeClusters
    file "seqs" from seqdb

    output:
    file "profile" into cascadeProfile

    script:
    QUERY = "seqs"
    TARGET = "seqs"
    RESULTS = "clusters"
    OUTDB = "profile"
    NCPUS = task.cpus
    template "mmseqs_result_to_profile.sh"
}


if ( params.enrich_profile ) {

    /*
     * Enrich the sequences by searching a database.
     */
    process enrichProfile {
        label "mmseqs"

        input:
        file "profile" from cascadeProfile
        file "enrich_seqs" from enrichdb

        output:
        file "enrich_matches" into enrichedSearchResults 

        script:
        PROFILE = "profile"
        TARGET = "enrich_seqs"
        OUTDB = "enrich_matches"
        NCPUS = task.cpus
        template "mmseqs_search_profile_strict.sh"
    }


    /*
     * Convert search results into an enriched profile.
     */
    process createEnrichedProfile {
        label "mmseqs"

        input:
        file "input_profile" from cascadeProfile
        file "enrich_seqs" from enrichdb
        file "enrich_matches" from enrichedSearchResults

        output:
        file "enriched_profile" into enrichedProfile

        script:
        QUERY = "input_profile"
        TARGET = "enrich_seqs"
        RESULTS = "enrich_matches"
        OUTDB = "enriched_profile"
        NCPUS = task.cpus
        template "mmseqs_result_to_profile.sh"
    }

    enrichedProfile.set { profile4Clu }
} else {
    cascadeProfile.set { profile4Clu }
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
    file "profile" into profileClusters
    file "profile_matches.tsv" into profileClustersResults

    script:
    INDB = "input_profile"
    OUTDB = "profile"
    NCPUS = task.cpus
    template "mmseqs_cluster_profile.sh"
}


/*
 * Merge the clustering results into single db.
 */
process mergeClusters {
    label "mmseqs"
    publishDir "${params.outdir}/clusters"

    input:
    file "seq" from seqdb
    file "cascade" from cascadeClusters
    file "profile" from profileClusters

    output:
    file "merged" into mergedClusters

    """
    mkdir -p merged

    mmseqs mergeclusters \
      seq/db \
      merged/db \
      cascade/db \
      profile/db
    """
}


/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractProfileClusterStats {
    label 'mmseqs'
    publishDir "${params.outdir}/clusters"

    input:
    file "profile" from mergedClusters
    file "seq" from seqdb

    output:
    file "profile.tsv" into profileClustersTSV
    file "profile_rep.fasta" into profileClustersRepFasta
    file "profile_stats.tsv" into profileClustersStats

    script:
    SEQS = "seq"
    CLUSTERS = "profile"
    NCPUS = task.cpus
    template "mmseqs_cluster_stats.sh"
}



if ( params.enrich_msa ) {
    /*
     * Create a profile database from the final clusters.
     */
    process createProfileCluProfile {
        label "mmseqs"

        input:
        file "clusters" from mergedClusters
        file "seqs" from seqdb

        output:
        file "profile" into profileClustersProfile

        script:
        QUERY = "seqs"
        TARGET = "seqs"
        RESULTS = "clusters"
        OUTDB = "profile"
	NCPUS = task.cpus
        template "mmseqs_result_to_profile.sh"
    }

    /*
     * Search the enrichment database to enrich msas.
     */
    process enrichMSAs {
        label "mmseqs"
        publishDir "msas"

        input:
        file "profile" from profileClustersProfile
        file "enrich_seqs" from enrichdb

        output:
        file "enrich_matches" into enrichMsaSearchResults 

        script:
        PROFILE = "profile"
        TARGET = "enrich_seqs"
        OUTDB = "enrich_matches"
        NCPUS = task.cpus
        template "mmseqs_search_profile_relaxed.sh"
    }


    process enrichMSAsResults {
        label "mmseqs"
        publishDir "msas"

        input:
        file "profile" from profileClustersProfile
        file "enrich_seqs" from enrichdb
        file "matches" from enrichMsaSearchResults

        output:
        file "enrich_matches.tsv" into enrichedMatches

        """
        mmseqs convertalis \
          profile/db \
          enrich_seqs/db \
          matches/db \
          "enrich_matches.tsv" \
          --threads ${task.cpus} \
          --format-mode 0 \
          --format-output "query,target,evalue,qcov,tcov,gapopen,pident,nident,mismatch,raw,bits,qstart,qend,tstart,tend,qlen,tlen,alnlen"

        sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' enrich_matches.tsv
        """
    }

    process mmseqsEnrichedMSA {
        label 'mmseqs'
        publishDir "msas"
        
        input:
        file "clusters" from enrichMsaSearchResults
        file "seq" from seqdb
        file "target" from enrichdb
        
        output:
        file "enriched" into enrichedMsa
        
        """
        mkdir -p enriched
        mmseqs result2msa \
          seq/db \
          target/db \
          clusters/db \
          enriched/db \
          --compress \
          --threads ${task.cpus}
        """
    }
}

if ( params.mafft ) {

    process clusterSeqdb {
        label "mmseqs"
    
        input:
        file "clusters" from mergedClusters
        file "seq" from seqdb
    
        output:
        set "cluster_seqs_*{,.index}" into splitClusters mode flatten
    
        """
        mkdir -p cluster_seqs
        mmseqs createseqfiledb "seq/db" "clusters/db" "cluster_seqs/db"
    
        NCLUSTERS=\$(wc -l < "cluster_seqs/db.index")
        NSPLITS=\$(( (\${NCLUSTERS} + 50000 + 1) / 50000 ))
    
        mmseqs splitdb "cluster_seqs/db" "cluster_seqs" --split \${NSPLITS}
        """
     
    } 
    
    
    process mafftMSA {
        label "mmseqs_mafft"
    
        input:
        set file(db), file(db_index) from splitClusters
            .map { [it.getSimpleName(), it] }
            .groupTuple(by: 0, size: 2)
            .map { bn, files -> files }
    
        output:
        set file("msas_${db.getName()}"), file("msas_${db_index.getName()}") into splitMuscleMSAs
    
        """
        mmseqs apply "${db}" "msas_${db.getName()}" -- mafft --retree 2 --maxiterate 2 --amino - 
        """
    }

    process combineSplitMSAs {
        label "mmseqs"
        publishDir "msas"
    
        input:
        file "*" from msas4CombineSplitMSAs
            .flatten()
            .collect()
    
        output:
        file "mafft" into combinedMSAs
    
        """
        mkdir -p mafft
        mmseqs concatdbs \$(find . -type f -not -name "*.index") mafft/db --preserve-keys
        """
    }

    msa4EstimateTrees = combinedMSAs
} else {
    process mmseqsMSA {
        label 'mmseqs'
        publishDir "msas"
        
        input:
        file "seq" from seqdb
        file "clusters" from mergedClusters
        
        output:
        file "mmseqs" into msa
        
        """
        mkdir -p mmseqs
        mmseqs result2msa \
          seq/db \
          seq/db \
          clusters/db \
          mmseqs/db \
          --compress \
          --threads ${task.cpus}

        cp seq/db_h mmseqs/db_header.ffdata
        cp seq/db_h.index mmseqs/db_header.ffindex
        cp seq/db mmseqs/db_sequence.ffdata
        cp seq/db.index mmseqs/db_sequence.ffindex
        """
    }
    msa4EstimateTrees = msa
}


if ( params.hhdata ) {

    process addCS219 {
        label "hhsuite"

        cpus 16
    
        input:
        file "msa" from msa
        file "hhdata" from hhdata
    
        output:
        file "hhmmseqs" into hhcs219
    
        """
        # Copy the directory because they all have to be kept together
        # and the alternative is to have some weird mutable folder thing happening.
        mkdir -p hhmmseqs
        cp -r msa/* hhmmseqs

        mpirun -np ${task.cpus} cstranslate_mpi \
            -i hhmmseqs/db \
            -o hhmmseqs/db_cs219 \
            -x 0.3 \
            -c 4 \
            -b \
            -I ca3m \
            -A hhdata/cs219.lib \
            -D hhdata/context_data.lib

        sort -k3 -n hhmmseqs/db_cs219.ffindex | cut -f1 > hhmmseqs/sorting.dat


        ffindex_order hhmmseqs/sorting.dat \
            hhmmseqs/db_ca3m.ff{data,index} \
            db_ca3m_ordered.ff{data,index}

        mv db_ca3m_ordered.ffdata hhmmseqs/db_ca3m.ffdata
        mv db_ca3m_ordered.ffindex hhmmseqs/db_ca3m.ffindex


        ffindex_order hhmmseqs/sorting.dat \
            hhmmseqs/db_cs219.ff{data,index} \
            db_cs219_ordered.ff{data,index}

        mv db_cs219_ordered.ffdata hhmmseqs/db_cs219.ffdata
        mv db_cs219_ordered.ffindex hhmmseqs/db_cs219.ffindex
        """
    }

    process addHHM {
        label "hhsuite"
        publishDir "msas"
        cpus 16

        input:
        file "hhcs219" from hhcs219

        output:
        file "hhmmseqs" into hhmsas

        """
        mkdir -p hhmmseqs
        ln -s \$(pwd)/hhcs219/* \$(pwd)/hhmmseqs

        a3m_database_extract \
            -i hhmmseqs/db_ca3m \
            -o db_a3m \
            -d hhmmseqs/db_sequence \
            -q hhmmseqs/db_header

        mpirun -np ${task.cpus} ffindex_apply_mpi \
            db_a3m.ff{data,index} \
            -i hhmmseqs/db_hhm.ffindex \
            -d hhmmseqs/db_hhm.ffdata \
            -- \
            hhmake -i stdin -o stdout -v 0

        ffindex_order hhmmseqs/sorting.dat \
            hhmmseqs/db_hhm.ff{data,index} \
            db_hhm_ordered.ff{data,index}

        mv db_hhm_ordered.ffdata hhmmseqs/db_hhm.ffdata
        mv db_hhm_ordered.ffindex hhmmseqs/db_hhm.ffindex

        rm db_a3m.ff{data,index}
        """
    }
}

    /*
if ( params.trees ) {
     * Estimate trees using the MSAs
    process estimateTrees {
        label "fasttree"
        publishDir "${params.outdir}/msas/trees"
        tag { msa.baseName }

        input:
        file msa from msas4EstimateTrees

        output:
        file "${msa.baseName}.nwk" into indivTrees

        """
        FastTree -fastest -quiet < "${msa}" > "${msa.baseName}.nwk"
        """
    }
}
     */

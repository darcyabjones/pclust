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
params.enrich_db = false
params.enrich_seqs = false
params.mafft = false
params.hhdata = false


if ( params.hhdata ) {
    hhdata = Channel
        .fromPath( params.hhdata, type: 'dir', checkIfExists: true, glob: false )
        .first()
} else {

    process getHHData {
        label "hhsuite"
        label "small_task"

        output:
        file "hhdata" into hhdata

        script:
        """
        cp -r \${HHLIB}/data hhdata
        """
    }
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
        """
        mkdir -p "seqdb"
        mmseqs createdb "seqs.fasta" "seqdb/db" --max-seq-len 14000

	# mkdir -p tmp
        # mmseqs createindex "seqdb/db" tmp --threads "${task.cpus}"
        # rm -rf -- tmp
        """
    }

} else {
    log.info "Please provide either sequences or the seqdb"
    exit 1
}


if ( params.enrich_db ) {

    enrichdb = Channel
        .fromPath( params.enrich_db, type: 'dir', checkIfExists: true, glob: false )
        .first()

} else if ( params.enrich_seqs ) {

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
        """
        mkdir -p "enrich_db"
        mmseqs createdb "seqs.fasta" "enrich_db/db" --max-seq-len 14000

        # mkdir -p tmp
        # mmseqs createindex "enrich_db/db" tmp --threads "${task.cpus}"
        # rm -rf -- tmp
        """
    }

} else {
    enrichdb = seqdb
}


/*
 * Perform the first pass clustering using basic mmseqs workflow.
 */
process clusterCascade {
    label 'mmseqs'
    label "big_task"

    input:
    file "seq" from seqdb

    output:
    file "cascade" into cascadeClusters

    script:
    INDB = "seq"
    OUTDB = "cascade"
    NCPUS = task.cpus

    """
    mkdir -p "${OUTDB}"
    mkdir -p tmp

    # mmseqs touchdb "${INDB}/db" --threads "${task.cpus}"

    mmseqs cluster \
      "${INDB}/db" \
      "${OUTDB}/db" \
      tmp \
      --threads "${NCPUS}" \
      --min-seq-id 0.0 \
      -c 0.8 \
      --cov-mode 0 \
      --cluster-steps 3 \
      -s 6.5 \
      --cluster-mode 0 \
      --db-load-mode 0

    rm -rf -- tmp
    """
}


/*
 * Extract information about each cluster.
 * E.G. cluster members, representative sequences, alignment statistics.
 */
process extractCascadeClusterStats {
    label 'mmseqs'
    label "big_task"

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


/*
 * Convert the cascade clusters into PSSMs
 */
process createProfile {
    label "mmseqs"
    label "big_task"
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


if ( params.enrich_db || params.enrich_seqs ) {

    /*
     * Enrich the sequences by searching a database.
     */
    process enrichProfile {
        label "mmseqs"
        label "big_task"

        input:
        file "profile" from cascadeProfile
        file "enrich_seqs" from enrichdb
	file "seqs" from seqdb

        output:
        file "enrich_matches" into enrichedSearchResults 

        script:
        """
	mkdir -p "tmp"
	mkdir -p "enrich_matches"

	mmseqs search \
	  "profile/db" \
	  "enrich_seqs/db" \
	  "enrich_matches/db" \
	  tmp \
	  --threads "${task.cpus}" \
	  --max-seqs 300 \
	  -e 0.00001 \
	  --e-profile 0.01 \
	  --start-sens 4.0 \
	  --sens-steps 3 \
	  -s 7.5 \
	  --rescore-mode 1 \
          --db-load-mode 0 \
          --split 0

	rm -rf -- tmp
        """
    }


    /*
     * Convert search results into an enriched profile.
     */
    process createEnrichedProfile {
        label "mmseqs"
        label "big_task"

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

    enrichedProfile.into { profile4CluSearch; profile4Clu }
} else {
    cascadeProfile.into { profile4CluSearch; profile4Clu }
}


/*
 * Perform the third clustering pass using sequence profiles.
 * This gets clusters down to about 10%-20% identity.
 */
process clusterProfileSearch {
    label 'mmseqs'
    label "big_task"

    publishDir "${params.outdir}/clusters"

    input:
    file "input_profile" from profile4CluSearch

    output:
    file "profile_matches" into profileClusterSearchResults
    file "profile_matches.tsv" into profileClusterSearchResultsTSV

    script:
    """
    mkdir -p tmp
    mkdir profile_matches
    
    mmseqs search \
      "input_profile/db" \
      "input_profile/db_consensus" \
      "profile_matches/db" \
      "tmp" \
      --threads "${task.cpus}" \
      --max-seqs 100 \
      -c 0.8 \
      --cov-mode 0 \
      --start-sens 4 \
      --sens-steps 3 \
      -s 7.5 \
      -e 0.00001 \
      --e-profile 0.01 \
      --db-load-mode 0 \
      --split 0 \
      --add-self-matches

    mmseqs convertalis \
      "input_profile/db" \
      "input_profile/db_consensus" \
      "profile_matches/db" \
      "profile_matches.tsv" \
      --threads "${task.cpus}" \
      --format-mode 0 \
      --format-output "query,target,evalue,qcov,tcov,gapopen,pident,nident,mismatch,raw,bits,qstart,qend,tstart,tend,qlen,tlen,alnlen"
    
    sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' "profile_matches.tsv"

    rm -rf -- tmp
    """
}


process clusterProfile {

    label 'mmseqs'
    label "big_task"

    publishDir "${params.outdir}/clusters"

    input:
    file "input_profile" from profile4Clu
    file "profile_matches" from profileClusterSearchResults

    output:
    file "profile" into profileClusters

    script:
    """
    mkdir -p "profile"
    mmseqs clust \
      "input_profile/db" \
      "profile_matches/db" \
      "profile/db" \
      --threads "${task.cpus}" \
      --cluster-mode 0
    """

}


/*
 * Merge the clustering results into single db.
 */
process mergeClusters {
    label "mmseqs"
    label "big_task"

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
    label "big_task"

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


/*
 */
process clusterSeqdb {
    label "mmseqs"
    label "big_task"

    input:
    file "clusters" from mergedClusters
    file "seq" from seqdb

    output:
    file "cluster_seqs_*" into splitClusters mode flatten

    script:
    """
    mkdir -p cluster_seqs
    mmseqs createseqfiledb "seq/db" "clusters/db" "cluster_seqs/db"

    TARGET_CLUSTER_SIZE=10000
    NCLUSTERS=\$(wc -l < "cluster_seqs/db.index")
    NSPLITS=\$(( (\${NCLUSTERS} + \${TARGET_CLUSTER_SIZE} + 1) / \${TARGET_CLUSTER_SIZE} ))

    mmseqs splitdb "cluster_seqs/db" "cluster_seqs" --split \${NSPLITS}
    """
}


process mafftMSA {
    label "mafft"
    label "big_task"

    input:
    set file(db), file(db_type), file(db_index)  from splitClusters
        .map { [it.getSimpleName(), it] }
        .groupTuple(by: 0, size: 3)
        .map { bn, files -> files }

    output:
    file "msas_${db.getName()}" into splitMafftMSAs

    script:
    """
    mkdir -p "msas_${db.getName()}"
    mpirun -np "${task.cpus}" mmseqs apply \
      "${db}" \
      "msas_${db.getName()}/db" \
      --threads 1 \
      -- \
      run_mafft.sh
    """
}

splitMafftMSAs.into {
    splitMafftMSAs4CombineSplitMSAs;
    splitMafftMSAs4EnrichMSA;
}


process combineSplitMSAs {
    label "mmseqs"
    label "small_task"

    publishDir "${params.outdir}/msas"

    input:
    file "*" from splitMafftMSAs4CombineSplitMSAs.collect()

    output:
    file "mafft" into combinedMSAs

    script:
    """
    for db in \$(find . -name "msas*")
    do
      if [ -e "mafft" ]
      then
        mkdir "reduce"
        mmseqs concatdbs "mafft/db" "\${db}/db" "reduce/db" --preserve-keys
        rm -rf -- "mafft"
        mv "reduce" "mafft"
      else
        cp -r -L "\${db}" "mafft"
      fi
    done
    """
}


/*
*/
process enrichMSA {

    label "mmseqs"
    label "big_task"

    input:
    file "msa" from splitMafftMSAs4EnrichMSA
    file "enrich" from enrichdb

    output:
    file "fas" into splitFas

    script:
    """
    mkdir profile
    mmseqs msa2profile msa/db profile/db --match-mode 1 --match-ratio 1 

    mkdir search tmp
    mmseqs search profile/db enrich/db search/db tmp -a 

    mkdir fas
    mmseqs result2msa profile/db enrich/db search/db fas/db

    mv fas/db fas/db.ffdata
    mv fas/db.dbtype fas/db.ffdata.dbtype
    mv fas/db.index fas/db.ffindex

    rm -rf -- profile search tmp
    """
}


/*
 */
process fasToHHDB {
    label "hhsuite"
    label "big_task"

    input:
    file "fas" from splitFas
    file "hhdata" from hhdata

    output:
    file "hhdb" into splitHHDB

    script:
    """
    mkdir -p hhdb
    cp fas/db.ffdata hhdb/db_fasta.ffdata
    cp fas/db.ffdata.dbtype hhdb/db_fasta.ffdata.dbtype
    cp fas/db.ffindex hhdb/db_fasta.ffindex

    mpirun -np "${task.cpus}" ffindex_apply_mpi \
        hhdb/db_fasta.ff{data,index} \
        -i hhdb/db_a3m.ffindex \
        -d hhdb/db_a3m.ffdata \
        -- \
        run_fas_to_a3m.sh

    mpirun -np "${task.cpus}" cstranslate_mpi \
        -i hhdb/db_a3m \
        -o hhdb/db_cs219 \
        -x 0.3 \
        -c 4 \
        -b \
        -I a3m \
        -A hhdata/cs219.lib \
        -D hhdata/context_data.lib

    mpirun -np "${task.cpus}" ffindex_apply_mpi \
        hhdb/db_a3m.ff{data,index} \
        -i hhdb/db_hhm.ffindex \
        -d hhdb/db_hhm.ffdata \
        -- \
        hhmake -i stdin -o stdout -v 0
    """
}


process combineSplitHHsuiteDbs {

    label "ffdb"
    label "small_task"

    publishDir "${params.outdir}"

    input:
    file "split_hhdata_*" from splitHHDB.collect()

    output:
    file "hhdata" into HHDB

    script:
    """
    mkdir unsorted hhdata

    ffdb combine \
      -d unsorted/db_fasta.ffdata \
      -i unsorted/db_fasta.ffindex \
      split_hhdata_*/db_fasta.ff{data,index}

    ffdb combine \
      -d unsorted/db_a3m.ffdata \
      -i unsorted/db_a3m.ffindex \
      split_hhdata_*/db_a3m.ff{data,index}

    ffdb combine \
      -d unsorted/db_cs219.ffdata \
      -i unsorted/db_cs219.ffindex \
      split_hhdata_*/db_cs219.ff{data,index}

    ffdb combine \
      -d unsorted/db_hhm.ffdata \
      -i unsorted/db_hhm.ffindex \
      split_hhdata_*/db_hhm.ff{data,index}

    sort -k3 -n unsorted/db_cs219.ffindex | cut -f1 > sorting.dat

    ffindex_order sorting.dat unsorted/db_cs219.ff{data,index} hhdata/db_cs219.ff{data,index}
    ffindex_order sorting.dat unsorted/db_a3m.ff{data,index} hhdata/db_a3m.ff{data,index}
    ffindex_order sorting.dat unsorted/db_fasta.ff{data,index} hhdata/db_fasta.ff{data,index}
    ffindex_order sorting.dat unsorted/db_hhm.ff{data,index} hhdata/db_hhm.ff{data,index}

    # rm -rf -- unsorted
    """
}

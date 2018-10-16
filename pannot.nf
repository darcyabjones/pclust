#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.seqdb = "$baseDir/clusters/sequence{,.dbtype,.index,.lookup,_h,_h.index}"

sequenceDB = Channel.fromPath( params.seqdb ).collect().into {
    sequenceDB1;
    sequenceDB2;
    sequenceDB3;
    sequenceDB4
}

process getUniqueSequences {
    container "soedinglab/mmseqs2"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
        "sequence_h", "sequence_h.index" from sequenceDB6

    output:
    file "unique.tsv" into uniqueProteinsTsv
    file "unique.fasta" into uniqueProteins

    """
    mkdir -p tmp
    mmseqs linclust sequence unique tmp --min-seq-id 1.0 -c 1.0

    mmseqs createtsv \
      sequence \
      sequence \
      unique \
      unique.tsv

    mmseqs result2repseq \
      sequence \
      unique \
      unique_rep

    mmseqs result2flat \
      sequence \
      sequence \
      unique_rep \
      unique_rep.fasta \
      --use-fasta-header

    mv unique_rep.fasta unique.fasta
    """
}

uniqueProteins.splitFasta(by: 2000).into {
    uniqueProteins1;
    uniqueProteins2;
    uniqueProteins3;
    uniqueProteins4;
    uniqueProteins5
    }

process effectorp {
    container "pclust/sperschneider"

    input:
    file fasta from uniqueProteins1

    output:
    file "table.tsv" into effectorpChunkedResults

    """
    EffectorP.py -s -i "${fasta}" > table.tsv
    """
}

effectorpResults = effectorpChunkedResults
    .collectFile(name: "annotations/effectorp.tsv")



/*
// Cluster a third time with HMMs.

process hmms {
    container "pclust/hhblits_mmseqs2"
    
    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB3
    set "profile_clusters", "profile_clusters.index" from profileClustDB1

    output:
    file "profile_clusters_hmms" into profileClustHMMs
    
    """
    mmseqs result2msa sequence sequence profile_clusters profile_clusters_msa --compress
    ln -s sequence_h profile_clusters_msa_header.ffdata
    ln -s sequence_h.index profile_clusters_msa_header.ffindex
    ln -s sequence profile_clusters_msa_sequence.ffdata
    ln -s sequence.index profile_clusters_msa_sequence.ffindex
    mpirun -np 2 cstranslate_mpi -i profile_clusters_msa -o profile_clusters_hmms -A /opt/hh-suite/data/cs219.lib -D /opt/hh-suite/data/context_data.lib -x 0.3 -c 4 -I ca3m
    """
}

*/

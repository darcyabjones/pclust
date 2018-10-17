#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.seqdb = "$baseDir/clusters/sequence"
seqdb = file(params.seqdb)
seqdbtype = file("${params.seqdb}.dbtype")
seqdbindex = file("${params.seqdb}.index")
seqdblookup = file("${params.seqdb}.lookup")
seqdbh = file("${params.seqdb}_h")
seqdbhindex = file("${params.seqdb}_h.index")

sequenceDB = Channel.create()
sequenceDB << [seqdb, seqdbtype, seqdbindex, seqdblookup, seqdbh, seqdbhindex]
sequenceDB.close()

sequenceDB.into {
    sequenceDB1;
    sequenceDB2;
    sequenceDB3;
    sequenceDB4
}

process getUniqueSequences {
    container "soedinglab/mmseqs2"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
        "sequence_h", "sequence_h.index" from sequenceDB1

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

uniqueProteins.splitFasta(by: 500).into {
    uniqueProteins1;
    uniqueProteins2;
    uniqueProteins3;
    uniqueProteins4;
    uniqueProteins5;
    uniqueProteins6;
    uniqueProteins7
    }

process effectorp {
    container "pclust/sperschneider"

    input:
    file fasta from uniqueProteins1

    output:
    file "${fasta}.tsv" into effectorpChunkedResults

    """
    EffectorP.py -s -i "${fasta}" -o "${fasta}.tsv"
    """
}

process gatherEffectorp {
    publishDir "annotations"

    input:
    file "tables" from effectorpChunkedResults.collect()

    output:
    file "effectorp.tsv" into effectorpResults

    """
    echo "seqid\teffector\tprobability" > effectorp.tsv
    cat tables* | grep -v "#" >> effectorp.tsv
    """
}

process signalp3 {
    container "pclust/signalp3"

    input:
    file fasta from uniqueProteins2

    output:
    file "${fasta}.tsv" into signalp3ChunkedResults

    """
    signalp -type euk -method "hmm" -short "${fasta}" > "${fasta}.tsv"
    """
}

process gatherSignalp3 {
    publishDir "annotations"

    input:
    file "tables" from signalp3ChunkedResults.collect()

    output:
    file "signalp3.tsv" into signalp3Results

    """
    echo "seqid\tsecreted\tcmax\tpos\tpos_decision\tsprob\tsprob_decision" > signalp3.tsv
    cat tables* | grep -v "#" | sed "s/ \\+/\t/g" >> signalp3.tsv
    """
}

process signalp4 {
    container "pclust/signalp4"

    input:
    file fasta from uniqueProteins3

    output:
    file "${fasta}.tsv" into signalp4ChunkedResults
    file "${fasta}_mature.fasta" into signalp4MatureProteins

    """
    signalp -t euk -f short -m "${fasta}_mature.fasta" "${fasta}" > "${fasta}.tsv"
    """
}

process gatherSignalp4 {
    publishDir "annotations"

    input:
    file "tables" from signalp4ChunkedResults.collect()

    output:
    file "signalp4.tsv" into signalp4Results

    """
    echo "seqid\tcmax\tcmax_pos\tymax\tymax_pos\tsmax\tsmax_pos\tsmean\td\tsecreted\tdmaxcut\tnetworks" > signalp4.tsv
    cat tables* | grep -v "#" | sed "s/ \\+/\t/g" >> signalp4.tsv
    """
}

process tmhmm {
    container "pclust/tmhmm"

    input:
    file fasta from uniqueProteins4

    output:
    file "${fasta}.tsv" into tmhmmChunkedResults

    """
    tmhmm -short -d < "${fasta}" > "${fasta}.tsv"
    """
}

process gatherTmhmm {
    publishDir "annotations"

    input:
    file "tables" from tmhmmChunkedResults.collect()

    output:
    file "tmhmm.tsv" into tmhmmResults

    """
    echo "seqid\tlen\texpaa\tfirst60\tpredhel\ttopology" > tmhmm.tsv

    cat tables* \
    | sed "s/len=\\|ExpAA=\\|First60=\\|PredHel=\\|Topology=//g" \
    | sed "s/ \\+/\t/g" \
    >> tmhmm.tsv
    """
}

process targetp {
    container "pclust/targetp"

    input:
    file fasta from uniqueProteins5

    output:
    file "${fasta}.tsv" into targetpChunkedResults

    """
    targetp -c -N < ${fasta} | tail -n+9 | head -n-2 > "${fasta}.tsv"
    """
}

process gatherTargetp {
    publishDir "annotations"

    input:
    file "tables" from targetpChunkedResults.collect()

    output:
    file "targetp.tsv" into targetpResults

    """
    echo "seqid\tlen\tmtp\tsp\tother\tloc\trc\ttplen" > targetp.tsv
    cat tables* | sed 's/ \\+/\t/g' >> targetp.tsv
    """
}

process targetpPlant {
    container "pclust/targetp"

    input:
    file fasta from uniqueProteins6

    output:
    file "${fasta}.tsv" into targetpPlantChunkedResults

    """
    targetp -c -P < ${fasta} | tail -n+9 | head -n-2 > "${fasta}.tsv"
    """
}

process gatherTargetpPlant {
    publishDir "annotations"

    input:
    file "tables" from targetpPlantChunkedResults.collect()

    output:
    file "targetp_plant.tsv" into targetpPlantResults

    """
    echo "seqid\tlen\tctp\tmtp\tsp\tother\tloc\trc\ttplen" > targetp_plant.tsv
    cat tables* | sed 's/ \\+/\t/g' >> targetp_plant.tsv
    """
}

process phobius {
    container "pclust/phobius"

    input:
    file fasta from uniqueProteins7

    output:
    file "${fasta}.tsv" into phobiusChunkedResults

    """
    sed 's/\\*\$//g' ${fasta} | phobius -short | tail -n+2 > "${fasta}.tsv"
    """
}

process gatherPhobius {
    publishDir "annotations"

    input:
    file "tables" from phobiusChunkedResults.collect()

    output:
    file "phobius.tsv" into phobiusResults

    """
    echo "seqid\ttm\tsp\tprediction" > phobius.tsv
    cat tables* | sed 's/ \\+/\t/g' >> phobius.tsv
    """
}

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

#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    =================================
    pclust/pannot
    =================================

    Usage:

    abaaab

    Mandatory Arguments:
      --seqs               description

    Options:
      --non-existant          description

    Outputs:

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}


// The unique proteins to run predictions on.
params.seqs = "$baseDir/sequences/dedup.fasta"
seqs = Channel.fromPath( params.seqs )


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
    mmseqs clusthash sequence/db result --threads ${task.cpus} --min-seq-id 1.0

    mkdir -p dedup
    mmseqs clust sequence/db result dedup/db --threads ${task.cpus}

    mmseqs createtsv sequence/db sequence/db dedup/db dedup.tsv --threads ${task.cpus}
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
    mmseqs result2repseq sequence/db cluster/db dedup --threads ${task.cpus}
    mmseqs result2flat sequence/db sequence/db dedup dedup.fasta --use-fasta-header
    """
}

/* Split the channel for reuse.
 * Note that targetp will be split separately into smaller bits because it's
 * temperamental.
 */
seqs.tap { seqs4Targetp }
    .splitFasta(by: 500)
    .into {
        seqs4Effectorp;
        seqs4Signalp3HMM;
        seqs4Signalp3NN;
        seqs4Signalp4;
        seqs4Tmhmm;
        seqs4Phobius;
        seqs4Apoplastp;
        seqs4LocalizerPlant;
    }


/*
 * Run effectorp on each sequence.
 */
process effectorp {
    label "sperschneider"

    input:
    file fasta from seqs4Effectorp

    output:
    file "${fasta}.tsv" into effectorpChunkedResults

    """
    EffectorP.py -s -i "${fasta}" -o "${fasta}.tsv"
    """
}


/*
 * Combine chunks of effectorp results into file.
 */
process gatherEffectorp {
    label "posix"
    publishDir "annotations"

    input:
    file "tables" from effectorpChunkedResults.collect()

    output:
    file "effectorp.tsv" into effectorpResults

    """
    echo "seqid\teffector\tprobability" > effectorp.tsv
    cat tables* \
    | grep -v "#" \
    | awk -F'\t' 'OFS="\t" { sub(/[[:space:]].*/, "", \$1); print \$1, \$2, \$3}' \
    >> effectorp.tsv
    """
}


/*
 * Run signalp3 HMM predictions for each sequence.
 * This is known to be more sensitive for detecting oomycete effectors.
 * See doi: 10.3389/fpls.2015.01168
 */
process signalp3HMM {
    label "signalp3"

    input:
    file fasta from seqs4Signalp3HMM

    output:
    file "${fasta}.tsv" into signalp3HMMChunkedResults

    """
    signalp -type euk -method "hmm" -short "${fasta}" > "${fasta}.tsv"
    """
}


/*
 * Combine signalp3 results into file.
 */
process gatherSignalp3HMM {
    label "posix"
    publishDir "annotations"

    input:
    file "tables" from signalp3HMMChunkedResults.collect()

    output:
    file "signalp3_hmm.tsv" into signalp3HMMResults

    """
    echo "seqid\tsecreted\tcmax\tpos\tpos_decision\tsprob\tsprob_decision" > signalp3_hmm.tsv
    cat tables* | grep -v "#" | sed "s/ \\+/\t/g" | sed "s/\t$//g" >> signalp3_hmm.tsv
    """
}


/*
 * Run signalp3 neural net predictions for each sequence.
 * This is known to be more sensitive for detecting fungal effectors.
 * See doi: 10.3389/fpls.2015.01168
process signalp3NN {
    label "signalp3"

    input:
    file fasta from seqs4Signalp3NN

    output:
    file "${fasta}.tsv" into signalp3NNChunkedResults

    """
    signalp -type euk -method "nn" -short "${fasta}" > "${fasta}.tsv"
    """
}
 */


/*
 * Combine signalp3 results into file.
process gatherSignalp3NN {
    publishDir "annotations"

    input:
    file "tables" from signalp3NNChunkedResults.collect()

    output:
    file "signalp3_nn.tsv" into signalp3NNResults

    """
    echo "seqid\tsecreted\tcmax\tpos\tpos_decision\tsprob\tsprob_decision" > signalp3_nn.tsv
    cat tables* | grep -v "#" | sed "s/ \\+/\t/g" >> signalp3_nn.tsv
    """
}
 */


/*
 * Run signalp4 for all sequences. Also get mature proteins for later use.
 * If there were no secreted proteins, don't emit fasta.
 */
process signalp4 {
    label "signalp4"

    input:
    file fasta from seqs4Signalp4

    output:
    file "${fasta}.tsv" into signalp4ChunkedResults
    file "${fasta}_mature.fasta" optional true into matureProteins

    """
    signalp -t euk -f short -m "${fasta}_mature.fasta" "${fasta}" > "${fasta}.tsv"

    if [ ! -s "${fasta}_mature.fasta" ]; then
      rm -f "${fasta}_mature.fasta"
    fi
    """
}


/*
 * Combine signalp4 chunks into file.
 */
process gatherSignalp4 {
    label "posix"
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


/*
 * Run tmhmm transmembrane domain prediction.
 */
process tmhmm {
    label "tmhmm"

    input:
    file fasta from seqs4Tmhmm

    output:
    file "${fasta}.tsv" into tmhmmChunkedResults

    """
    tmhmm -short -d < "${fasta}" > "${fasta}.tsv"
    """
}


/*
 * Collect results into file
 */
process gatherTmhmm {
    label "posix"
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


/*
 * Because targetp is a bit finicky, process it in smaller chunks.
 */
seqs4Targetp
    .splitFasta(by: 100)
    .into {
        seqs4Targetp1;
        seqs4Targetp2;
    }


/*
 * Run targetp using non-plant networks for chunks.
 */
process targetp {
    label "targetp"

    input:
    file fasta from seqs4Targetp1

    output:
    file "${fasta}.tsv" into targetpChunkedResults

    """
    targetp -c -N < ${fasta} | tail -n+9 | head -n-2 > "${fasta}.tsv"
    """
}


/*
 * Collect targetp chunks into file
 */
process gatherTargetp {
    label "posix"
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


/*
 * Run targetp using plant networks for chunks.
 */
process targetpPlant {
    label "targetp"

    input:
    file fasta from seqs4Targetp2

    output:
    file "${fasta}.tsv" into targetpPlantChunkedResults

    """
    targetp -c -P < ${fasta} | tail -n+9 | head -n-2 > "${fasta}.tsv"
    """
}


/*
 * Collect targetp chunks into file.
 */
process gatherTargetpPlant {
    label "posix"
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


/*
 * Run phobius predictions for chunks.
 * Phobius has comparable sensitivity to signalp nn models and also runs tm prediction.
 */
process phobius {
    label "phobius"

    input:
    file fasta from seqs4Phobius

    output:
    file "${fasta}.tsv" into phobiusChunkedResults

    """
    sed 's/\\*\$//g' "${fasta}" | phobius -short | tail -n+2 > "${fasta}.tsv"
    """
}


/*
 * Collect phobius results into file.
 */
process gatherPhobius {
    label "posix"
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
 * Run apoplastp for chunks
 */
process apoplastp {
    label "sperschneider"

    input:
    file fasta from seqs4Apoplastp

    output:
    file "${fasta}.tsv" into apoplastpChunkedResults

    """
    ApoplastP.py -s -i "${fasta}" -o "${fasta}.tsv"
    """
}


/*
 * Collect chunked apoplastp results into file.
 */
process gatherApoplastp {
    label "posix"
    publishDir "annotations"

    input:
    file "tables" from apoplastpChunkedResults.collect()

    output:
    file "apoplastp.tsv" into apoplastpResults

    """
    echo "seqid\tprediction\tprobability" > apoplastp.tsv
    cat tables* \
    | grep -v "#" \
    | awk -F'\t' 'OFS="\t" { sub(/[[:space:]].*/, "", \$1); print \$1, \$2, \$3}' \
    >> apoplastp.tsv
    """
}


/*
 * Run localizer in "effector" mode using mature peptides from signalp
 */
process localizerEffector {
    label "sperschneider"

    input:
    file fasta from matureProteins

    output:
    file "${fasta}.tsv" into localizerEffectorChunkedResults

    """
    LOCALIZER.py -e -M -i "${fasta}" -o results
    grep -v "^#" results/Results.txt \
    | tail -n+2 \
    | sed '/^\\s*\$/d' \
    | awk -F'\t' 'OFS="\t" { sub(/[[:space:]].*/, "", \$1); print \$1, \$2, \$3, \$4}' \
    > "${fasta}.tsv"
    """
}


/*
 * Collect localizer results into a file.
 */
process gatherLocalizerEffector {
    label "posix"
    publishDir "annotations"

    input:
    file "tables" from localizerEffectorChunkedResults.collect()

    output:
    file "localizer_effector.tsv" into localizerEffectorResults

    """
    echo "seqid\tchloroplast\tmitochondria\tnucleus" > localizer_effector.tsv
    cat tables* >> localizer_effector.tsv
    """
}


/*
 * Run localizer using plant mode.
 */
process localizerPlant {
    label "sperschneider"

    input:
    file fasta from seqs4LocalizerPlant

    output:
    file "${fasta}.tsv" into localizerPlantChunkedResults

    """
    LOCALIZER.py -p -i "${fasta}" -o results
    grep -v "^#" results/Results.txt \
    | tail -n+2 \
    | sed '/^\\s*\$/d' \
    | awk -F'\t' 'OFS="\t" { sub(/[[:space:]].*/, "", \$1); print \$1, \$2, \$3, \$4}' \
    > "${fasta}.tsv"
    """
}


/*
 * Collect localizer results into file.
 */
process gatherLocalizerPlant {
    label "posix"
    publishDir "annotations"

    input:
    file "tables" from localizerPlantChunkedResults.collect()

    output:
    file "localizer_plant.tsv" into localizerPlantResults

    """
    echo "seqid\tchloroplast\tmitochondria\tnucleus" > localizer_plant.tsv
    cat tables* >> localizer_plant.tsv
    """
}


/*
 * Combine annotation results into single file.
 */
process combineAnnotations {
    label "R"
    publishDir "annotations"

    input:
    file "apoplastp.tsv" from apoplastpResults
    file "effectorp.tsv" from effectorpResults
    file "localizer_effector.tsv" from localizerEffectorResults
    file "localizer_plant.tsv" from localizerPlantResults
    file "signalp3_hmm.tsv" from signalp3HMMResults
    file "signalp4.tsv" from signalp4Results
    file "targetp_non_plant.tsv" from targetpResults
    file "targetp_plant.tsv" from targetpPlantResults
    file "tmhmm.tsv" from tmhmmResults
    file "phobius.tsv" from phobiusResults

    output:
    file "combined.tsv" into combinedResults

    """
    join_annotations.R \
      apoplastp.tsv \
      effectorp.tsv \
      localizer_effector.tsv \
      localizer_plant.tsv \
      signalp3_hmm.tsv \
      dummy \
      signalp4.tsv \
      targetp_non_plant.tsv \
      targetp_plant.tsv \
      tmhmm.tsv \
      phobius.tsv \
    > combined.tsv
    """
}

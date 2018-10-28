#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.seqs = "$baseDir/dedup/dedup.fasta"
seqs = Channel.fromPath( params.seqs )

seqs.tap { seqs4Targetp }
    .splitFasta(by: 500)
    .into {
        seqs4Effectorp;
        seqs4Signalp3;
        seqs4Signalp4;
        seqs4Tmhmm;
        seqs4Phobius;
        seqs4Apoplastp;
        seqs4LocalizerPlant;
    }

process effectorp {
    container "pclust/sperschneider"

    input:
    file fasta from seqs4Effectorp

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
    cat tables* \
    | grep -v "#" \
    | awk -F'\t' 'OFS="\t" { sub(/[[:space:]].*/, "", \$1); print \$1, \$2, \$3}' \
    >> effectorp.tsv
    """
}

process signalp3 {
    container "pclust/signalp3"

    input:
    file fasta from seqs4Signalp3

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
    file fasta from seqs4Signalp4

    output:
    file "${fasta}.tsv" into signalp4ChunkedResults
    file "${fasta}_mature.fasta" into matureProteins

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
    file fasta from seqs4Tmhmm

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

seqs4Targetp
    .splitFasta(by: 100)
    .into {
        seqs4Targetp1;
        seqs4Targetp2;
    }

process targetp {
    container "pclust/targetp"

    input:
    file fasta from seqs4Targetp1

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
    file fasta from seqs4Targetp2

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
    file fasta from seqs4Phobius

    output:
    file "${fasta}.tsv" into phobiusChunkedResults

    """
    sed 's/\\*\$//g' "${fasta}" | phobius -short | tail -n+2 > "${fasta}.tsv"
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

process apoplastp {
    container "pclust/sperschneider"

    input:
    file fasta from seqs4Apoplastp

    output:
    file "${fasta}.tsv" into apoplastpChunkedResults

    """
    ApoplastP.py -s -i "${fasta}" -o "${fasta}.tsv"
    """
}

process gatherApoplastp {
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

process localizerEffector {
    container "pclust/sperschneider"

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

process gatherLocalizerEffector {
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

process localizerPlant {
    container "pclust/sperschneider"

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

process gatherLocalizerPlant {
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

#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.msas = "$baseDir/msas/muscle_refine/*.faa"
msas = Channel.fromPath(params.msas)

process createHmms {
    container "pclust/hhblits_mmseqs2"
    publishDir "hmms"

    input:
    file msa from msas

    output:
    file "${msa.baseName}.hhm" into hmms

    """
    hhmake \
      -i "${msa}" \
      -o "${msa.baseName}.hhm" \
      -id 90 \
      -M first 
    """
}

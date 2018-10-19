#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.msas = "$baseDir/msas/muscle_refine/*.faa"
msas = Channel.fromPath(params.msas)

params.uniref90 = none
// Download uniref90 if not provided
// Download dssp rsync -avz rsync://rsync.cmbi.ru.nl/dssp/ ./dssp
// Download pdb
// Download hhblits databases

// Add secondary structure info using psipred

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

//

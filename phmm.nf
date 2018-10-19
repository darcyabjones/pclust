#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.msas = "$baseDir/msas/muscle_refine/*.faa"
msas = Channel.fromPath(params.msas)

params.dssp = false
params.pdb = false
params.hhpdb = false
params.hhscop = false
params.hhuniclust = false
params.hhpfam = false

// Download dssp rsync -avz rsync://rsync.cmbi.ru.nl/dssp/ ./dssp
// Download pdb
// Download hhblits databases
//          - uniclust 30 http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz
//          - scop http://wwwuser.gwdg.de/%7Ecompbiol/data/hhsuite/databases/hhsuite_dbs/scop90_01Mar17.tgz
//          - pdb http://wwwuser.gwdg.de/%7Ecompbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_181017.tar.gz
//          - pfam http://wwwuser.gwdg.de/%7Ecompbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_31.0.tgz


process enrichMsas {
    container "pclust/hhblits"

    input:
    file msa from msas
    file db from hhuniclust

    output:
    file "${msa.baseName}.a3m" into hmms

    """
    hhblits \
      -i "${msa}" \
      -oa3m "${msa.baseName}.a3m" \
      -id 90 \
      -cov 60 \
      -M 50 \
      -mact 0.3 \
      -n 2 \
      -d "${db}"
    """
}

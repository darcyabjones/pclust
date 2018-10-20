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
params.hhuniref = false
params.hhpfam = false


if ( !params.pdb ) {
    process downloadPDB {
        label "download"
        storeDir "databases"

        output:
        file "pdb" into pdbDatabase

        """
        rsync \
          -rlpt \
          -v \
          -z \
          --delete \
          --port=33444 \
          rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ \
          ./pdb
        """
    }
} else if ( file(params.pdb).exists() ) {
    dsspDatabase = Channel.fromPath( params.pdb )
} else {
    exit 1, "You specified a pdb database, but it doesn't exist."
}


if ( !params.dssp ) {
    process downloadDSSP {
        label "download"
        storeDir "databases"

        output:
        file "dssp" into dsspDatabase

        """
        rsync -avz rsync://rsync.cmbi.ru.nl/dssp/ ./dssp
        """
    }
} else if ( file(params.dssp).exists() ) {
    dsspDatabase = Channel.fromPath( params.dssp )
} else {
    exit 1, "You specified a dssp database, but it doesn't exist."
}


if ( !params.hhuniref ) {
    process downloadHHUniref {
        label "download"
        storeDir "databases"

        output:
        file "hhuniref" into hhunirefDatabase

        """
        mkdir -p hhuniref
        cd hhuniref
        wget -c http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz
        tar -zxf uniclust30_2018_08_hhsuite.tar.gz
        mv uniclust30_2018_08/* ./
        rm -rf -- uniref
        """
    }
} else if ( file(params.hhuniref).exists() ) {
    hhunirefDatabase = Channel.fromPath( params.hhuniref )
} else {
    exit 1, "You specified a hhblits formatted uniref database, but it doesn't exist."
}


if ( !params.hhpdb ) {
    process downloadHHPDB {
        label "download"
        storeDir "databases"

        output:
        file "hhpdb" into hhpdbDatabase

        """
        mkdir -p hhpdb
        cd hhpdb
        wget -c http://wwwuser.gwdg.de/%7Ecompbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_181017.tar.gz
        tar -zxf pdb70_from_mmcif_181017.tar.gz
        """
    }
} else if ( file(params.hhpdb).exists() ) {
    hhpdbDatabase = Channel.fromPath( params.hhpdb )
} else {
    exit 1, "You specified a hhblits formatted pdb database, but it doesn't exist."
}


if ( !params.hhscop ) {
    process downloadHHSCOP {
        label "download"
        storeDir "databases"

        output:
        file "hhscop" into hhscopDatabase

        """
        mkdir -p hhscop
        cd hhscop
        wget -c http://wwwuser.gwdg.de/%7Ecompbiol/data/hhsuite/databases/hhsuite_dbs/scop90_01Mar17.tgz
        tar -zxf scop90_01Mar17.tgz
        """
    }
} else if ( file(params.hhscop).exists() ) {
    hhscopDatabase = Channel.fromPath( params.hhscop )
} else {
    exit 1, "You specified a hhblits formatted scop database, but it doesn't exist."
}


if ( !params.hhpfam ) {
    process downloadHHPFAM {
        label "download"
        storeDir "databases"

        output:
        file "hhpfam" into hhpfamDatabase

        """
        mkdir -p hhpfam
        cd hhpfam
        wget -c http://wwwuser.gwdg.de/%7Ecompbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_31.0.tgz
        tar -zxf pfamA_31.0.tgz
        """
    }
} else if ( file(params.hhpfam).exists() ) {
    hhpfamDatabase = Channel.fromPath( params.hhpfam )
} else {
    exit 1, "You specified a hhblits formatted pfam database, but it doesn't exist."
}


process enrichMsas {
    container "pclust/hhblits"

    input:
    file msa from msas
    file db from hhunirefDatabase

    output:
    file "${msa.baseName}.a3m" into hmms

    """
    hhblits \
      -i "${msa}" \
      -o "${msa}.hhr" \
      -oa3m "${msa.baseName}.a3m" \
      -ohhm "${msa.baseName}.hhm" \
      -id 90 \
      -cov 60 \
      -M 50 \
      -mact 0.3 \
      -n 2 \
      -d "${db}/uniclust30_2018_08"
    """
}


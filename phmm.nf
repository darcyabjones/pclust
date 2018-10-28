#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    =================================
    pclust/phmm
    =================================

    Usage:

    abaaab

    Mandatory Arguments:
      --genomes               description

    Options:
      --non-existant          description

    Outputs:

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

params.msas = "$baseDir/msas/muscle_refine/*.faa"
msas = Channel.fromPath(params.msas)

params.dssp = false
params.pdb = false
params.hhpdb = false
params.hhscop = false
params.hhuniref = false
params.hhpfam = false

/*
Get the databases ready.
Downloads them if they don't exist already.
*/

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
        rm -rf -- uniclust30_2018_08
        """
    }
} else if ( file(params.hhuniref).exists() ) {
    hhunirefDatabase = Channel.fromPath( params.hhuniref )
} else {
    exit 1, "You specified a hhblits formatted uniref database, but it doesn't exist."
}

hhunirefDatabase.into { hhunirefDatabase4Enrich; hhunirefDatabase4Search }


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


/*
 * Enrich the MSAs using hmms from uniref.
 */
process enrichMsas {
    label "hhblits"
    publishDir "msas/enriched"

    input:
    file msa from msas
    file "db" from hhunirefDatabase4Enrich

    output:
    file "${msa.baseName}.hhr" into enrichedResults
    file "${msa.baseName}.a3m" into enrichedMsas

    """
    NSEQS=\$(grep -c "^>" "${msa}")
    if [ "\${NSEQS}" -eq "1" ]; then
      MOPT="first"
    else
      MOPT="50"
    fi

    hhblits \
      -i "${msa}" \
      -o "${msa.baseName}.hhr" \
      -oa3m "${msa.baseName}.a3m" \
      -ohhm "${msa.baseName}.hhm" \
      -atab "${msa.baseName}.tsv" \
      -id 90 \
      -cov 60 \
      -M \${MOPT} \
      -mact 0.4 \
      -cpu 4 \
      -d "db/uniclust30_2018_08"
    """
}

enrichedMsas.into {
    msas4Database;
    msas4Search;
    msas4Uniref;
    msas4Pfam;
    msas4Scop;
    msas4Pdb
}


/*
 * Combine the enriched MSAs into a database that can be searched with hhsearch
 */
process createHmmDatabase {
    label "hhblits"
    publishDir "hhsuite"

    input:
    file "*.a3m" from msas4Database.collect()

    output:
    file "clusterdb" into clusterDb

    """
    mkdir -p clusterdb
    hhsuitedb.py \
      --ia3m *.a3m \
      -o clusterdb/db \
      --cpu 4
    """
}

process searchClusters {
    label "hhblits"
    publishDir "hhsuite/clusters"

    input:
    file msa from msas4Search
    file "db" from clusterDb

    output:
    file "${msa.baseName}_result.hhr" into clusterResults
    file "${msa.baseName}_result.a3m" into clusterMsas

    """
    hhsearch \
      -i "${msa}" \
      -o "${msa.baseName}_result.hhr" \
      -oa3m "${msa.baseName}_result.a3m" \
      -atab "${msa.baseName}_result.tsv" \
      -M a2m \
      -all \
      -mact 0.4 \
      -cpu 4 \
      -d "db/db"
    """
}


process searchUniref {
    label "hhblits"
    publishDir "hhsuite/uniref"

    input:
    file msa from msas4Uniref
    file "db" from hhunirefDatabase4Search

    output:
    file "${msa.baseName}_result.hhr" into unirefResults
    file "${msa.baseName}_result.a3m" into unirefMsas

    """
    hhblits \
      -i "${msa}" \
      -o "${msa.baseName}_result.hhr" \
      -oa3m "${msa.baseName}_result.a3m" \
      -atab "${msa.baseName}_result.tsv" \
      -n 1 \
      -M a2m \
      -all \
      -mact 0.4 \
      -cpu 4 \
      -d db/uniclust30_2018_08
    """
}


process searchPfam {
    label "hhblits"
    publishDir "hhsuite/pfam"

    input:
    file msa from msas4Pfam
    file "db" from hhpfamDatabase

    output:
    file "${msa.baseName}_result.hhr" into pfamResults
    file "${msa.baseName}_result.a3m" into pfamMsas

    """
    hhsearch \
      -i "${msa}" \
      -o "${msa.baseName}_result.hhr" \
      -oa3m "${msa.baseName}_result.a3m" \
      -atab "${msa.baseName}_result.tsv" \
      -M a2m \
      -all \
      -mact 0.4 \
      -cpu 4 \
      -d db/pfam
    """
}

process searchScop {
    label "hhblits"
    publishDir "hhsuite/scop"

    input:
    file msa from msas4Scop
    file "db" from hhscopDatabase

    output:
    file "${msa.baseName}_result.hhr" into scopResults
    file "${msa.baseName}_result.a3m" into scopMsas

    """
    hhsearch \
      -i "${msa}" \
      -o "${msa.baseName}_result.hhr" \
      -oa3m "${msa.baseName}_result.a3m" \
      -atab "${msa.baseName}_result.tsv" \
      -M a2m \
      -all \
      -mact 0.4 \
      -cpu 4 \
      -d "db/scop90"
    """
}

process searchPdb {
    label "hhblits"
    publishDir "hhsuite/pdb"

    input:
    file msa from msas4Pdb
    file "db" from hhpdbDatabase

    output:
    file "${msa.baseName}_result.hhr" into pdbResults
    file "${msa.baseName}_result.a3m" into pdbMsas

    """
    hhsearch \
      -i "${msa}" \
      -o "${msa.baseName}_result.hhr" \
      -oa3m "${msa.baseName}_result.a3m" \
      -atab "${msa.baseName}_result.tsv" \
      -M a2m \
      -all \
      -mact 0.4 \
      -cpu 4 \
      -d "db/pdb70"
    """
}

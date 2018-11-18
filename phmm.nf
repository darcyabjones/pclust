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
      --msas <glob pattern> Multiple sequence alignements to analyse.

    Options:
      --enrich              Enrich input MSAs using uniref. 
      --nopdb               Don't search for matches in PDB
      --nopfam              Don't search for matches in PFam
      --noscop              Don't search for matches in SCOP
      --nouniref            Don't search for matches in uniref

    Database download options:
      --hhpdb              Use this directory containing hhsuite PDB database.
      --hhpfam             Use this directory containing hhsuite PFam database.
      --hhscop             Use this directory containing hhsuite SCOP database.
      --hhuniref           Use this directory containing hhsuite Uniref database.

      Any of these options will prevent download of db and use path instead.

    Outputs:

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

params.msas = false

params.enrich = false
params.nopdb = false
params.nopfam = false
params.noscop = false
params.nouniref = false

//params.dssp = false
//params.pdb = false

params.hhpdb = false
params.hhscop = false
params.hhuniref = false
params.hhpfam = false


if ( !params.msas ) {
    log.info "Hey I need some MSAs to use please."
    exit 1
}

msas = Channel
    .fromPath(params.msas)
    .map { [it.baseName, it] }

/*
Get the databases ready.
Downloads them if they don't exist already.
*/

/*
if ( !params.pdb ) {
    process downloadPDB {
        label "download"
        storeDir "${params.outdir}/databases"

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
        storeDir "${params.outdir}/databases"

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
*/

if ( !params.hhuniref ) {
    process downloadHHUniref {
        label "download"
        storeDir "${params.outdir}/databases"

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
        storeDir "${params.outdir}/databases"

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
        storeDir "${params.outdir}/databases"

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
        storeDir "${params.outdir}/databases"

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


if ( params.enrich ) {
    /*
     * Enrich the MSAs using hmms from uniref.
     */
    process enrichMsas {
        label "hhblits"
        publishDir "${params.outdir}/enriched_msas"
        tag { label }
    
        input:
        set val(label), file("input.faa") from msas
        file "db" from hhunirefDatabase4Enrich
    
        output:
        set val(label), file("${label}.hhr") into enrichedResults
        set val(label), file("${label}.a3m") into enrichedMsas
    
        """
        NSEQS=\$(grep -c "^>" "input.faa")
        if [ "\${NSEQS}" -eq "1" ]; then
          MOPT="first"
        else
          MOPT="50"
        fi
    
        hhblits \
          -i "input.faa" \
          -o "${label}.hhr" \
          -oa3m "${label}.a3m" \
          -ohhm "${label}.hhm" \
          -atab "${label}.tsv" \
          -n 3 \
          -id 90 \
          -M \${MOPT} \
          -mact 0.4 \
          -maxmem 12.0 \
          -cpu ${task.cpus} \
          -d "db/uniclust30_2018_08"
        """
    }
} else {
    /*
     * Convert fasta MSAs to a3m format
     */
    process msa2A3m {
        label "hhblits"
        tag { label }

        input:
        set val(label), file("input.faa") from msas

        output:
        set val(label), file("${label}.a3m") into enrichedMsas

        """
        NSEQS=\$(grep -c "^>" "input.faa")
        if [ "\${NSEQS}" -eq "1" ]; then
          MOPT="first"
        else
          MOPT="50"
        fi

        reformat.pl \
          fas \
          a3m \
          "input.faa" \
          "${label}.a3m" \
          -M \${MOPT}
        """
    }
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
 * Combine the enriched MSAs into a database that can be searched with hhsearch.
 */
process createHmmDatabase {
    label "hhblits"
    publishDir "${params.outdir}/hhsuite"

    input:
    file "*.a3m" from msas4Database.map {l, f -> f} .collect()

    output:
    file "clusterdb" into clusterDb

    """
    mkdir -p clusterdb
    hhsuitedb.py \
      --ia3m *.a3m \
      -o clusterdb/db \
      --cpu ${task.cpus}
    """
}


/*
 * Search each cluster hmm against all other clusters.
 */
process searchClusters {
    label "hhblits"
    publishDir "${params.outdir}/hhsuite/clusters"
    tag { label }

    input:
    set val(label), file("input.a3m") from msas4Search
    file "db" from clusterDb

    output:
    set val(label), file("${label}.hhr") into clusterResults
    set val(label), file("${label}.a3m") into clusterMsas

    """
    hhsearch \
      -i "input.a3m" \
      -o "${label}.hhr" \
      -oa3m "${label}.a3m" \
      -atab "${label}.tsv" \
      -M a2m \
      -all \
      -mact 0.4 \
      -cpu ${task.cpus} \
      -d "db/db"
    """
}


if ( !params.nouniref ) {
    /*
     * Search for shorter matches in the uniref database.
     */
    process searchUniref {
        label "hhblits"
        publishDir "${params.outdir}/hhsuite/uniref"
        tag { label }
    
        input:
        set val(label), file("input.a3m") from msas4Uniref
        file "db" from hhunirefDatabase4Search
    
        output:
        set val(label), file("${label}.hhr") into unirefResults
        set val(label), file("${label}.a3m") into unirefMsas
    
        """
        hhblits \
          -i "input.a3m" \
          -o "${label}.hhr" \
          -oa3m "${label}.a3m" \
          -atab "${label}.tsv" \
          -n 1 \
          -M a2m \
          -all \
          -mact 0.4 \
          -cpu ${task.cpus} \
          -d db/uniclust30_2018_08
        """
    }
}


if ( !params.nopfam ) {
    /*
     * Search for pfam domains.
     */
    process searchPfam {
        label "hhblits"
        publishDir "${params.outdir}/hhsuite/pfam"
        tag { label }
    
        input:
        set val(label), file("input.a3m") from msas4Pfam
        file "db" from hhpfamDatabase
    
        output:
        set val(label), file("${label}.hhr") into pfamResults
        set val(label), file("${label}.a3m") into pfamMsas
    
        """
        hhsearch \
          -i "input.a3m" \
          -o "${label}.hhr" \
          -oa3m "${label}.a3m" \
          -atab "${label}.tsv" \
          -M a2m \
          -all \
          -mact 0.4 \
          -cpu ${task.cpus} \
          -d db/pfam
        """
    }
}


if ( !params.noscop ) {
    /*
     * Search for SCOP matches.
     */
    process searchScop {
        label "hhblits"
        publishDir "${params.outdir}/hhsuite/scop"
        tag { label }
    
        input:
        set val(label), file("input.a3m") from msas4Scop
        file "db" from hhscopDatabase
    
        output:
        set val(label), file("${label}.hhr") into scopResults
        set val(label), file("${label}.a3m") into scopMsas
    
        """
        hhsearch \
          -i "input.a3m" \
          -o "${label}.hhr" \
          -oa3m "${label}.a3m" \
          -atab "${label}.tsv" \
          -M a2m \
          -all \
          -mact 0.4 \
          -cpu ${task.cpus} \
          -d "db/scop90"
        """
    }
}


if ( !params.nopdb ) {
    /*
     * Search for PDB matches.
     */
    process searchPdb {
        label "hhblits"
        publishDir "${params.outdir}/hhsuite/pdb"
        tag { label }
    
        input:
        set val(label), file("input.a3m") from msas4Pdb
        file "db" from hhpdbDatabase
    
        output:
        set val(label), file("${label}.hhr") into pdbResults
        set val(label), file("${label}.a3m") into pdbMsas
    
        """
        hhsearch \
          -i "input.a3m" \
          -o "${label}.hhr" \
          -oa3m "${label}.a3m" \
          -atab "${label}.tsv" \
          -M a2m \
          -all \
          -mact 0.4 \
          -cpu ${task.cpus} \
          -d "db/pdb70"
        """
    }
}

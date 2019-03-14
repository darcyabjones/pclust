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
      --clusters <glob pattern>

    Options:
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

params.clusters = false

params.nopdb = false
params.nopfam = false
params.noscop = false
params.nouniref = false

params.hhpdb = false
params.hhscop = false
params.hhuniref = false
params.hhpfam = false


if ( params.clusters ) {
    clusters = Channel
        .fromPath(params.clusters, type: 'dir', checkIfExists: true, glob: false)
        .first()
} else {
    log.info "Hey I need some MSAs to use please."
    exit 1
}

/*
 * Get the databases ready.
 * Downloads them if they don't exist already.
 */


if ( !params.nouniref && params.hhuniref ) {
    hhunirefDatabase = Channel
        .fromPath( params.hhuniref, type: 'dir', checkIfExists: true, glob: false)
        .first()
} else if ( !params.nouniref ) {
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
}


if ( !params.nopdb && params.hhpdb ) {
    hhpdbDatabase = Channel
        .fromPath( params.hhpdb, type: 'dir', checkIfExists: true, glob: false )
        .first()
} else if ( !params.nopdb ) {
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
}

if ( !params.noscop && params.hhscop ) {
    hhscopDatabase = Channel
        .fromPath( params.hhscop, type: 'dir', checkIfExists: true, glob: false )
        .first()
} else if ( !params.noscop ) {
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
}


if ( !params.nopfam && params.hhpfam ) {
    hhpfamDatabase = Channel
        .fromPath( params.hhpfam, type: 'dir', checkIfExists: true, glob: false )
        .first()
} else if ( !params.nopfam ) {
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
}


/*
 * Run the database searches.
 */

process decompressDB {
    label "hhsuite"

    input:
    file "clusters" from clusters

    output:
    file "decompressed" into decompressed

    """
    mkdir -p decompressed

    a3m_database_extract \
        -i clusters/db_ca3m \
        -o decompressed/db_a3m \
        -d clusters/db_sequence \
        -q clusters/db_header
    """
}


process splitDB {
    label "python"

    input:
    file "clusters" from decompressed

    output:
    // NB we rely on the data, index order of this glob, throughout.
    // Changing it will cause errors.
    file "subset_*.ff{data,index}" into splitClusters mode flatten

    """
    ffdb.py split -n 10000 -b "subset_{index}.{ext}" clusters/db_a3m.ff{data,index}
    """
}


splitClusters
    .map { [it.getSimpleName(), it] }
    .groupTuple(by: 0, size: 2)
    .map { bn, files -> [bn, files[0], files[1]] }
    .into {
        clusters4Search;
        clusters4Uniref;
        clusters4Pfam;
        clusters4Scop;
        clusters4Pdb;
    }



/*
 * Search each cluster hmm against all other clusters.
 */
process searchClusters {
    label "hhsuite"
    tag { name }
    cpus = 16

    input:
    set val(name), file("subset.ffdata"), file("subset.ffindex") from clusters4Search
    file "clusters" from clusters

    output:
    set val("clusters"), file("${name}.ffdata"), file("${name}.ffindex") into clusterSearchResults

    script:
    CPU_PER_TASK = 4
    NTASKS = task.cpus.intdiv(CPU_PER_TASK)
    INPUT = "subset"
    OUTPUT = name
    DB  = "clusters"

    template "hhblits_mpi_sensitive.sh"
}


if ( !params.nouniref ) {
    /*
     * Search for shorter matches in the uniref database.
     */
    process searchUniref {
        label "hhblits"
        tag { name }
    
        input:
        set val(name), file("subset.ffdata"), file("subset.ffindex") from clusters4Uniref
        file "uniref" from hhunirefDatabase
    
        output:
        set val("uniref"), file("${name}.ffdata"), file("${name}.ffindex") into unirefSearchResults
    
        script:
        CPU_PER_TASK = 2
        NTASKS = task.cpus.intdiv(CPU_PER_TASK)
        INPUT = "subset"
        OUTPUT = name
        DB  = "uniref"
        NITER = 1
    
        template "hhblits_mpi.sh"
    }
} else {
    unirefSearchResults = Channel.empty()
}


if ( !params.nopfam ) {
    /*
     * Search for pfam domains.
     */
    process searchPfam {
        label "hhblits"
        tag { name }
    
        input:
        set val(name), file("subset.ffdata"), file("subset.ffindex") from clusters4Pfam
        file "pfam" from hhpfamDatabase
    
        output:
        set val("pfam"), file("${name}.ffdata"), file("${name}.ffindex") into pfamSearchResults
    
        script:
        CPU_PER_TASK = 2
        NTASKS = task.cpus.intdiv(CPU_PER_TASK)
        INPUT = "subset"
        OUTPUT = name
        DB  = "pfam"
    
        template "hhblits_mpi_sensitive.sh"
    }
} else {
    pfamSearchResults = Channel.empty()
}


if ( !params.noscop ) {
    /*
     * Search for SCOP matches.
     */
    process searchScop {
        label "hhblits"
        tag { name }
    
        input:
        set val(name), file("subset.ffdata"), file("subset.ffindex") from clusters4Scop
        file "scop" from hhscopDatabase
    
        output:
        set val("scop"), file("${name}.ffdata"), file("${name}.ffindex") into scopSearchResults
    
        script:
        CPU_PER_TASK = 2
        NTASKS = task.cpus.intdiv(CPU_PER_TASK)
        INPUT = "subset"
        OUTPUT = name
        DB  = "scop"
    
        template "hhblits_mpi_sensitive.sh"
    }
} else {
    scopSearchResults = Channel.empty()
}


if ( !params.nopdb ) {
    /*
     * Search for PDB matches.
     */
    process searchPdb {
        label "hhblits"
        tag { name }
    
        input:
        set val(name), file("subset.ffdata"), file("subset.ffindex") from clusters4Pdb
        file "pdb" from hhpdbDatabase
    
        output:
        set val("pdb"), file("${name}.ffdata"), file("${name}.ffindex") into pdbSearchResults
    
    
        script:
        CPU_PER_TASK = 2
        NTASKS = task.cpus.intdiv(CPU_PER_TASK)
        INPUT = "subset"
        OUTPUT = name
        DB  = "pdb"
    
        template "hhblits_mpi_sensitive.sh"
    }
} else {
    pdbSearchResults = Channel.empty()
}


/*
 * Collect the search results into a single database.
 */
process collectSearchResults {
    label "python"
    tag { db }

    input:
    set val(db), file("*") from clusterSearchResults
        .concat(unirefSearchResults, pfamSearchResults, scopSearchResults, pdbSearchResults)
        .flatMap { db, data, index -> [[db, data], [db, index]] }
        .groupTuple(by: 0) 

    output:
    set val(db), file("${db}_results.ffdata"), file("${db}_results.ffindex") into searchResults

    """
    ffdb.py combine -d "${db}_results.ffdata" -i "${db}_results.ffindex" *.{ffdata,ffindex} 
    """
}

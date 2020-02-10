#!/usr/bin/env nextflow


/*
 * Get fasta sequence databases from the cluster database and
 * and split the database into many parts to allow checkpointing.
 */
process split_db {

    label "mmseqs"
    label "small_task"

    input:
    val chunk_size
    val dbtype
    path "db"

    output:
    path "split_db_*", emit: chunks

    script: 
    """
    TARGET_CHUNK_SIZE="${chunk_size}"
    NENTRIES=\$(wc -l < "db/db.index")
    NSPLITS=\$(( (\${NENTRIES} + \${TARGET_CHUNK_SIZE} + 1) / \${TARGET_CHUNK_SIZE} ))

    cp -rL db db_tmp
    if [ "${dbtype}" = "protein" ]
    then
        awk 'BEGIN { printf("%c%c%c%c",0,0,0,0); exit; }' > db_tmp/db.dbtype
    elif [ "${dbtype}" = "msa" ]
    then
        awk 'BEGIN { printf("%c%c%c%c",11,0,0,0); exit; }' > db_tmp/db.dbtype
    elif [ "${dbtype}" = "nucleotide" ]
    then
        awk 'BEGIN { printf("%c%c%c%c",1,0,0,0); exit; }' > db_tmp/db.dbtype
    fi
    

    mmseqs splitdb "db_tmp/db" "tmp_split_db" --split \${NSPLITS}

    rm -rf -- db_tmp

    for f in tmp_split_db_*.index
    do
        BASENAME="\${f%.index}"
        DIRNAME="\${BASENAME#tmp_}"

        mkdir "\${DIRNAME}"
        mv "\${f}" "\${DIRNAME}/db.index"
        mv "\${BASENAME}" "\${DIRNAME}/db"
        mv "\${BASENAME}.dbtype" "\${DIRNAME}/db.dbtype"
    done
    """
}


/*
 * Collect split databases into a single database.
 */
process combine_split_dbs {

    label "mmseqs"
    label "small_task"

    input:
    val name
    path "split_db_*"

    output:
    path "${name}", emit: combined

    script:
    """
    for db in split_db_*
    do
      if [ -e "${name}" ]
      then
        mkdir "reduce"
        mmseqs concatdbs "${name}/db" "\${db}/db" "reduce/db" --preserve-keys
        rm -rf -- "${name}"
        mv "reduce" "${name}"
      else
        cp -r -L "\${db}" "${name}"
      fi
    done

    ORIG=\${PWD}
    cd "${name}"
    ln -s "db" "db.ffdata"
    ln -s "db.dbtype" "db.ffdata.dbtype"
    ln -s "db.index" "db.ffindex"
    cd "\${ORIG}"
    """
}


process check_user_db {

    label "mmseqs"
    label "small_task"

    input:
    path "user"

    output:
    path "checked", emit: checked

    script:
    """
    cp -r -L user checked

    ORIG=\${PWD}

    if [ ! -e checked/db ] && [ -e checked/db.ffdata ]
    then
        cd checked
        ln -s db.ffdata db
        cd "\${ORIG}"
    elif [ ! -e checked/db.ffdata ] && [ -e checked/db ]
    then
        cd checked
        ln -s db db.ffdata
        cd "\${ORIG}"
    elif [ -e checked/db ] && [ -e checked/db.ffdata ]
    then
        true
    else
        echo "A database was specified but there was no db or db.ffdata file!" 1>&2
        exit 1
    fi

    if [ ! -e checked/db.index ] && [ -e checked/db.ffindex ]
    then
        cd checked
        ln -s db.ffindex db.index
        cd "\${ORIG}"
    elif [ ! -e checked/db.ffindex ] && [ -e checked/db.index ]
    then
        cd checked
        ln -s db.index db.ffindex
        cd "\${ORIG}"
    elif [ -e checked/db ] && [ -e checked/db.ffdata ]
    then
        true
    else
        echo "A database was specified but there was no db.index or db.ffindex file!" 1>&2
        exit 1
    fi
    """
}


/*
 * Construct MSA from the fasta databases.
 */
process mafft {

    label "mafft"
    label "big_task"

    input:
    path "fastas"

    output:
    path "msas", emit: msas

    script:
    """
    mkdir -p "msas"
    mpirun -np "${task.cpus}" mmseqs apply \
      "fastas/db" \
      "msas/db" \
      --threads 1 \
      -- \
      run_mafft.sh
    """
}


/*
 * Construct an ML tree for each cluster.
 */
process fasttree {

    label "fasttree"
    label "big_task"

    input:
    path "msas"

    output:
    path "trees", emit: trees

    script:
    """
    OMP_NUM_THREADS=1

    mkdir -p "trees"
    mpirun -np "${task.cpus}" mmseqs apply \
      "msas/db" \
      "trees/db" \
      --threads 1 \
      -- \
      FastTree -quiet -nopr
    """
}


workflow split_user_db {

    get:
    chunk_size // Should be an integer value channel
    dbtype
    user

    main:
    checked = check_user_db(user)
    split_dbs = split_db(chunk_size, db_type, checked).chunks.flatten()

    emit:
    chunks
}

#!/usr/bin/env nextflow

process create_hhdb {

    label "hhsuite"
    label "big_task"

    input:
    path "msas"
    path "hhdata"

    output:
    path "hhdb", emit: hhdb

    script:
    """
    mkdir -p hhdb

    #for f in msas/*
    #do
    #    EXT="\$(basename \${f##*.})"
    #    cp -L "\${f}" "hhdb/db_fasta.\${EXT}"
    #done

    mpirun -np ${task.cpus} ffindex_apply_mpi \
        msas/db.ff{data,index} \
        -i hhdb/db_a3m.ffindex \
        -d hhdb/db_a3m.ffdata \
        -- \
        run_fas_to_a3m.sh

    mpirun -np "${task.cpus}" cstranslate_mpi \
        -i "hhdb/db_a3m" \
        -o "hhdb/db_cs219" \
        -x 0.3 \
        -c 4 \
        -b \
        -I a3m \
        -A "hhdata/cs219.lib" \
        -D "hhdata/context_data.lib"

    mpirun -np "${task.cpus}" ffindex_apply_mpi \
        hhdb/db_a3m.ff{data,index} \
        -i "hhdb/db_hhm.ffindex" \
        -d "hhdb/db_hhm.ffdata" \
        -- \
        hhmake -i stdin -o stdout -v 0
    """
}


process copy_hhdata {

    label "hhsuite"
    label "small_task"

    output:
    path "hhdata", emit: hhdata

    script:
    """
    cp -r \${HHLIB}/data hhdata
    """
}


process find_hhdbs_order {

    label "posix"
    label "small_task"

    input:
    path "split_hhdbs_*"

    output:
    path "order.txt", emit: order

    script:
    """
    cat split_hhdbs_*/db_cs219.ffindex \
    | sort -k3,3n \
    | cut -f1 \
    > order.txt
    """
}


process combine_sort_split_hhdbs {

    label "hhsuite"
    label "small_task"

    input:
    val kind
    path "split_hhdbs_*"
    path "order.txt"

    output:
    tuple path("db_${kind}.ffdata"),
          path("db_${kind}.ffindex"), emit: combined

    script:
    assert kind in ["cs219", "a3m", "hhm"]

    """
    ffdb combine \
        -d "db_tmp.ffdata" \
        -i "db_tmp.ffindex" \
        split_hhdbs_*/db_${kind}.ff{data,index}

    ffindex_order \
      "order.txt" \
      db_tmp.ff{data,index} \
      db_${kind}.ff{data,index}

    rm -f db_tmp.*
    """
}


process combine_sorted_hhdbs {

    label "posix"
    label "small_task"

    input:
    tuple path("db_cs219.ffdata"),
          path("db_cs219.ffindex")

    tuple path("db_a3m.ffdata"),
          path("db_a3m.ffindex")

    tuple path("db_hhm.ffdata"),
          path("db_hhm.ffindex")

    output:
    path "hhself", emit: db

    script:
    """
    mkdir hhself
    ln -sf "\${PWD}/"*.ff{data,index} "\${PWD}/hhself"
    """
}


process search_hmms {

    label "hhsuite"
    label "big_task"

    input:
    path "query"
    path "target"
    val nrealign

    output:
    path "matches"

    script:
    """
    mkdir "matches"
    mpirun -np "${task.cpus}" hhblits_mpi \
      -i "subdb/db_hhm" \
      -d "hhdb/db" \
      -o "results/db_hhr" \
      -n 1 \
      -cpu 1 \
      -v 0 \
      -e 0.001 \
      -E 0.001 \
      -z 0 \
      -Z "${nrealign}" \
      -b 0 \
      -B "${nrealign}" \
      -pre_evalue_thresh 10 \
      -min_prefilter_hits 10 \
      -realign_max "${nrealign}"
    """
}


process split_hhm_databases {

    label "ffdb"
    label "small_task"

    input:
    val chunk_size
    path "hhdb"

    output:
    path "split_db_*"

    script:
    """
    ffdb split \
        --size "${chunk_size}" \
        --basename "split_db_{index}/db.{ext}" \
        hhdb/db_hhm.ff{data,index}
    """
}


process combine_hhsuite_results {

    label "ffdb"
    label "small_task"

    input:
    path "split_results_*"

    output:
    path "combined"

    script:
    """
    mkdir -p "combined"

    ffdb combine \
      -d "combined/db_hhr.ffdata" \
      -i "combined/db_hhr.ffindex" \
      split_results_*/db_hhr.ff{data,index}
    """
}


// This is a way to deal with nf not liking to use the same process twice.
workflow combine_sort_split_cs219_hhdbs {

    get:
    split_hhdbs
    order

    main:
    combined = combine_sort_split_hhdbs("cs219", split_hhdbs.collect(), order)

    emit:
    combined
}


workflow combine_sort_split_a3m_hhdbs {

    get:
    split_hhdbs
    order

    main:
    combined = combine_sort_split_hhdbs("a3m", split_hhdbs.collect(), order)

    emit:
    combined
}


workflow combine_sort_split_hhm_hhdbs {

    get:
    split_hhdbs
    order

    main:
    combined = combine_sort_split_hhdbs("hhm", split_hhdbs.collect(), order)

    emit:
    combined
}


workflow create_and_sort_hhdb {


    get:
    msas // These should be split already
    hhdata

    main:
    split_hhdbs = create_hhdb(msas, hhdata)
    order = find_hhdbs_order(split_hhdbs.collect())
    hhdb_cs219 = combine_sort_split_cs219_hhdbs(split_hhdbs, order)
    hhdb_a3m = combine_sort_split_a3m_hhdbs(split_hhdbs, order)
    hhdb_hhm = combine_sort_split_hhm_hhdbs(split_hhdbs, order)
    hhdb = combine_sorted_hhdbs(hhdb_cs219, hhdb_a3m, hhdb_hhm)

    emit:
    hhdb
}

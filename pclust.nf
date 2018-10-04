#!/usr/bin/env nextflow

params.fastas = "$baseDir/data/*.faa"

fastas = Channel.fromPath( params.fastas )

// Need step to rename fastas to include filename.

process combineFasta {

    input:
    file "*.fa*" from fastas.collect()

    output:
    file "combined.faa" into combinedFasta

    """
    cat *.fa* > combined.faa
    """
}


sequenceDB.into { sequenceDB1; sequenceDB2; sequenceDB3; sequenceDB4; sequenceDB5; sequenceDB6; sequenceDB7; sequenceDB8; sequenceDB9; sequenceDB10; sequenceDB11 }

process mmseqsDB {
    container "soedinglab/mmseqs2"

    input:
    file fasta from combinedFasta

    output:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  into sequenceDB

    """
    mmseqs createdb "${fasta}" sequence --max-seq-len 14000
    """
}


process mmseqsLinclust {
    container "soedinglab/mmseqs2"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB1

    output:
    file "clusters", "clusters.index" into clustDB

    """
    mkdir -p tmp
    mmseqs linclust sequence clusters tmp
    """
}


process mmseqsExtractClusters {
    container: "soedinglab/mmseqs2"

    publishDir "clusters"

    input:
    set "clusters", "clusters.index" from clustDB
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB2

    output:
    file "clusters.tsv" into uniclust30TSV
    file "clusters_rep.fasta" into uniclust30Fasta

    """
    mmseqs createtsv \
      sequence \
      sequence \
      clusters \
      clusters.tsv

    mmseqs result2repseq \
      sequence \
      clusters \
      clusters_rep

    mmseqs result2flat \
      sequence \
      sequence \
      clusters_rep \
      clusters_rep.fasta \
      --use-fasta-header
    """
}


/*


sequenceDB.into { sequenceDB1; sequenceDB2; sequenceDB3; sequenceDB4; sequenceDB5; sequenceDB6; sequenceDB7; sequenceDB8; sequenceDB9; sequenceDB10; sequenceDB11 }


process mmseqsFragPref {
    container "soedinglab/mmseqs2"
    publishDir "test"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB1

    output:
    set "pref_frag", "pref_frag.index" into prefFrag

    """
    mmseqs prefilter \
      sequence \
      sequence \
      pref_frag \
      --max-seqs 4000 \
      --min-ungapped-score 100 \
      --comp-bias-corr 0 \
      -s 1
    """
}


prefFrag.into { prefFrag1; prefFrag2 }

process mmseqsFragDiag {
    container: "soedinglab/mmseqs2"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB2
    set "pref_frag", "pref_frag.index" from prefFrag1

    output:
    set "aln_frag", "aln_frag.index" into alnFrag

    """
    mmseqs rescorediagonal \
      sequence \
      sequence \
      pref_frag \
      aln_frag \
      --min-seq-id 0.9 \
      -c 0.95 \
      --cov-mode 1
    """
}


process mmseqsFragClust {
    container: "soedinglab/mmseqs2"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB3
    set "aln_frag", "aln_frag.index" from alnFrag

    output:
    set "clu_frag", "clu_frag.index" into cluFrag

    """
    mmseqs clust \
      sequence \
      aln_frag \
      clu_frag \
      --cluster-mode 2
    """
}

cluFrag.into { cluFrag1; cluFrag2; cluFrag3 }


process mmseqsFragSubDB {
    container: "soedinglab/mmseqs2"

    input:
    set "clu_frag", "clu_frag.index" from cluFrag1
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB4

    output:
    set "input_step_redundancy", "input_step_redundancy.dbtype", "input_step_redundancy.index" into fragDB 
    file "order_frag" into orderFrag

    """
    awk '{ print \$1 }' clu_frag.index > order_frag

    mmseqs createsubdb order_frag sequence input_step_redundancy
    """
}

fragDB.into { fragDB1; fragDB2; fragDB3; fragDB4 }


process mmseqsFilterFragHash {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step_redundancy", "input_step_redundancy.dbtype", "input_step_redundancy.index" from fragDB1
    
    output:
    set "aln_redundancy", "aln_redundancy.index" into alnRedundancy
    
    """
    mmseqs clusthash \
      input_step_redundancy \
      aln_redundancy \
      --min-seq-id 0.9
    """
}


process mmseqsFilterFragClust {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step_redundancy", "input_step_redundancy.dbtype", "input_step_redundancy.index" from fragDB2
    set "aln_redundancy", "aln_redundancy.index" from alnRedundancy
    
    output:
    set "clu_redundancy", "clu_redundancy.index" into cluRedundancy

    """
    mmseqs clust \
      input_step_redundancy \
      aln_redundancy \
      clu_redundancy \
      --cluster-mode 2
    """
}

cluRedundancy.into { cluRedundancy1; cluRedundancy2; cluRedundancy3; cluRedundancy4 }

process mmseqsFilterFragSubDB {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step_redundancy", "input_step_redundancy.dbtype", "input_step_redundancy.index" from fragDB3
    set "clu_redundancy", "clu_redundancy.index" from cluRedundancy1

    output:
    set "input_step0", "input_step0.dbtype", "input_step0.index" into step0DB
    file "order_redundancy" into orderRedundancy
    
    """
    awk '{ print \$1 }' clu_redundancy.index > order_redundancy
    
    mmseqs createsubdb \
      order_redundancy \
      input_step_redundancy \
      input_step0
    """
}

step0DB.into { step0DB1; step0DB2; step0DB3; step0DB4 }
orderRedundancy.into { orderRedundancy1; orderRedundancy2 }

process mmseqsCluster90SequenceFiltered {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step0", "input_step0.dbtype", "input_step0.index" from step0DB1
    file "order_redundancy" from orderRedundancy1
    set "pref_frag", "pref_frag.index" from prefFrag2

    output:
    set "pref_frag_filtered", "pref_frag_filtered.index" into prefFragFiltered

    """
    mmseqs createsubdb \
      order_redundancy \
      pref_frag \
      pref_frag_filtered
    """
}


process mmseqsCluster90FilterDB {
    container: "soedinglab/mmseqs2"

    input:
    set "pref_frag_filtered", "pref_frag_filtered.index" from prefFragFiltered
    file "order_redundancy" from orderRedundancy2

    output:
    set "pref_0", "pref_0.index" into pref0

    """
    mmseqs filterdb \
      pref_frag_filtered \
      pref_0 \
      --filter-file order_redundancy
    """
}


process mmseqsCluster90Align {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step0", "input_step0.dbtype", "input_step0.index" from step0DB2
    set "pref_0", "pref_0.index" from pref0

    output:
    set "aln_0", "aln_0.index" into aln0

    """
    mmseqs align \
      input_step0 \
      input_step0 \
      pref_0 \
      aln_0 \
      --max-seqs 100 \
      -c 0.9 \
      --alignment-mode 2 \
      --min-seq-id 0.9 \
      --comp-bias-corr 0 \
      -e 0.001 \
      --max-seq-len 32768 \
      --max-rejected 2147483647
    """
}


process mmseqsCluster90Clust {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step0", "input_step0.dbtype", "input_step0.index" from step0DB3
    set "aln_0", "aln_0.index" from aln0

    output:
    set "clu_0", "clu_0.index" into clu0

    """
    mmseqs clust \
      input_step0 \
      aln_0 \
      clu_0 \
      --cluster-mode 2
    """
}

clu0.into { clu0_1; clu0_2; clu0_3; clu0_4 }

process mmseqsCluster90SubDB {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step0", "input_step0.dbtype", "input_step0.index" from step0DB4
    set "clu_0", "clu_0.index" from clu0_1
    
    output:
    set "input_step1", "input_step1.dbtype", "input_step1.index" into step1DB
    file "order_0" into order0

    """
    awk '{ print \$1 }' clu_0.index > order_0

    mmseqs createsubdb order_0 input_step0 input_step1
    """
}

step1DB.into { step1DB1; step1DB2; step1DB3; step1DB4 }


// Step 1

process mmseqsCluster90BigPreFilterDB {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step1", "input_step1.dbtype", "input_step1.index" from step1DB1

    output:
    set "pref_1", "pref_1.index" into pref1

    """
    mmseqs prefilter \
      input_step1 \
      input_step1 \
      pref_1 \
      --max-seqs 100  \
      -c 0.9 \
      --comp-bias-corr 1 \
      -s 2
    """
}


process mmseqsCluster90BigAlign {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step1", "input_step1.dbtype", "input_step1.index" from step1DB2
    set "pref_1", "pref_1.index" from pref1

    output:
    set "aln_1", "aln_1.index" into aln1

    """
    mmseqs align \
      input_step1 \
      input_step1 \
      pref_1 \
      aln_1 \
      --max-seqs 100 \
      -c 0.8 \
      --alignment-mode 2 \
      --min-seq-id 0.9 \
      --comp-bias-corr 1 \
      -e 0.001 \
      --max-seq-len 32768 \
      --max-rejected 2147483647
    """
}


process mmseqsCluster90BigClust {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step1", "input_step1.dbtype", "input_step1.index" from step1DB3
    set "aln_1", "aln_1.index" from aln1

    output:
    set "clu_1", "clu_1.index" into clu1

    """
    mmseqs clust \
      input_step1 \
      aln_1 \
      clu_1 \
      --cluster-mode 0
    """
}

clu1.into { clu1_1; clu1_2; clu1_3; clu1_4 }


process mmseqsCluster90BigMerge {
    container: "soedinglab/mmseqs2"

    publishDir "uniclust90"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB6
    set "clu_frag", "clu_frag.index" from cluFrag2
    set "clu_redundancy", "clu_redundancy.index" from cluRedundancy2
    set "clu_0", "clu_0.index" from clu0_2
    set "clu_1", "clu_1.index" from clu1_1

    output:
    set "uniclust90", "uniclust90.index" into uniclust90

    """
    mmseqs mergeclusters \
      sequence \
      uniclust90 \
      clu_frag \
      clu_redundancy \
      clu_0 \
      clu_1
    """    
}


process mmseqsCluster90BigSubDB {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step1", "input_step1.dbtype", "input_step1.index" from step1DB4
    set "clu_1", "clu_1.index" from clu1_2

    output:
    file "order_1" into order1
    set "input_step2", "input_step2.dbtype", "input_step2.index" into step2DB

    """
    awk '{ print \$1 }' clu_1.index > order_1

    mmseqs createsubdb \
      order_1 \
      input_step1 \
      input_step2
    """
}

step2DB.into { step2DB1; step2DB2; step2DB3; step2DB4 }

process mmseqsClusterLowPrefilter {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step2", "input_step2.dbtype", "input_step2.index" from step2DB1

    output:
    set "pref_step2", "pref_step2.index" into pref2

    """
    mmseqs prefilter \
      input_step2 \
      input_step2 \
      pref_step2 \
      --max-seqs 300 \
      -c 0.8 \
      --comp-bias-corr 1 \
      -s 6 
    """
}


process mmseqsClusterLowAlign {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step2", "input_step2.dbtype", "input_step2.index" from step2DB2
    set "pref_step2", "pref_step2.index" from pref2

    output:
    set "aln_step2", "aln_step2.index" into aln2

    """
    mmseqs align \
      input_step2 \
      input_step2 \
      pref_step2 \
      aln_step2 \
      --max-seqs 300 \
      -c 0.8 \
      --alignment-mode 3 \
      --min-seq-id 0.3 \
      --comp-bias-corr 1 \
      -e 0.001 \
      --max-seq-len 32768 \
      --max-rejected 2147483647
    """
}

aln2.into { aln2_1; aln2_2 }

process mmseqsCluster50Filter {
    container: "soedinglab/mmseqs2"

    input:
    set "aln_step2", "aln_step2.index" from aln2_1

    output:
    set "aln_uniclust50", "aln_uniclust50.index" into alnUniclust50

    """
    mmseqs filterdb \
      aln_step2 \
      aln_uniclust50 \
      --filter-column 3 \
      --filter-regex '(0\\.[5-9][0-9]{2}|1\\.000)'
    """
}


process mmseqsCluster50Clust {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step2", "input_step2.dbtype", "input_step2.index" from step2DB3
    set "aln_uniclust50", "aln_uniclust50.index" from alnUniclust50

    output:
    set "clu_uniclust50", "clu_uniclust50.index" into cluUniclust50

    """
    mmseqs clust \
      input_step2 \
      aln_uniclust50 \
      clu_uniclust50 \
      --cluster-mode 0
    """
}


process mmseqsCluster50Merge {
    container: "soedinglab/mmseqs2"

    publishDir "uniclust50"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB7
    set "clu_frag", "clu_frag.index" from cluFrag2
    set "clu_redundancy", "clu_redundancy.index" from cluRedundancy3
    set "clu_0", "clu_0.index" from clu0_3
    set "clu_1", "clu_1.index" from clu1_2
    set "clu_uniclust50", "clu_uniclust50.index" from cluUniclust50

    output:
    set "uniclust50", "uniclust50.index" into uniclust50

    """
    mmseqs mergeclusters \
      sequence \
      uniclust50 \
      clu_frag \
      clu_redundancy \
      clu_0 \
      clu_1 \
      clu_uniclust50
    """
}


process mmseqsCluster30Clust {
    container: "soedinglab/mmseqs2"

    input:
    set "input_step2", "input_step2.dbtype", "input_step2.index" from step2DB4
    set "aln_step2", "aln_step2.index" from aln2_2
    
    output:
    set "clu_uniclust30", "clu_uniclust30.index" into cluUniclust30

    """
    mmseqs clust \
      input_step2 \
      aln_step2 \
      clu_uniclust30 \
      --cluster-mode 0
    """
}


process mmseqsCluster30Merge {
    container: "soedinglab/mmseqs2"

    publishDir "uniclust30"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB8
    set "clu_frag", "clu_frag.index" from cluFrag3
    set "clu_redundancy", "clu_redundancy.index" from cluRedundancy4
    set "clu_0", "clu_0.index" from clu0_4
    set "clu_1", "clu_1.index" from clu1_3
    set "clu_uniclust30", "clu_uniclust30.index" from cluUniclust30

    output:
    set "uniclust30", "uniclust30.index" into uniclust30

    """
    mmseqs mergeclusters \
      sequence \
      uniclust30 \
      clu_frag \
      clu_redundancy \
      clu_0 \
      clu_1 \
      clu_uniclust30
    """
}


process mmseqsCluster30TSV {
    container: "soedinglab/mmseqs2"

    publishDir "uniclust30"

    input:
    set "uniclust30", "uniclust30.index" from uniclust30
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB9

    output:
    file "uniclust30.tsv" into uniclust30TSV
    file "uniclust30_rep.fasta" into uniclust30Fasta

    """
    mmseqs createtsv \
      sequence \
      sequence \
      uniclust30 \
      uniclust30.tsv
    
    mmseqs result2repseq \
      sequence \
      uniclust30 \
      uniclust30_rep

    mmseqs result2flat \
      sequence \
      sequence \
      uniclust30_rep \
      uniclust30_rep.fasta \
      --use-fasta-header
    """
}

process mmseqsCluster50TSV {
    container: "soedinglab/mmseqs2"

    publishDir "uniclust50"

    input:
    set "uniclust50", "uniclust50.index" from uniclust30
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB10

    output:
    file "uniclust50.tsv" into uniclust50TSV

    """
    mmseqs createtsv \
      sequence \
      sequence \
      uniclust50 \
      uniclust50.tsv
    """
}

process mmseqsCluster90TSV {
    container: "soedinglab/mmseqs2"

    publishDir "uniclust90"

    input:
    set "uniclust90", "uniclust90.index" from uniclust90
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB11

    output:
    file "uniclust90.tsv" into uniclust90TSV

    """
    mmseqs createtsv \
      sequence \
      sequence \
      uniclust90 \
      uniclust90.tsv
    """
}
*/

//process mmseqsCluster30ToMSA {

//    publishDir ""
//}

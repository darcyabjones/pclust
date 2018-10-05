#!/usr/bin/env nextflow

params.fastas = "$baseDir/data/*.faa"

fastas = Channel.fromPath( params.fastas )


// Combine the fasta files in preparation for clustering
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


// Create the sequence database
// This will get reused a lot.
process createSequenceDB {
    container "soedinglab/mmseqs2"

    input:
    file fasta from combinedFasta

    output:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  into sequenceDB

    """
    mmseqs createdb "${fasta}" sequence --max-seq-len 14000
    """
}

// Duplicate the channel
sequenceDB.into { sequenceDB1; sequenceDB2; sequenceDB3; sequenceDB4;
                  sequenceDB5; sequenceDB6; sequenceDB7; sequenceDB8;
                  sequenceDB9; sequenceDB10; sequenceDB11 }


// Do the first clustering pass.
// This tends to get to about 50% identity.
process mmseqsClust {
    container "soedinglab/mmseqs2"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB1

    output:
    set "clusters", "clusters.index" into clustDB

    """
    mkdir -p tmp
    mmseqs cluster sequence clusters tmp -s 5 --cluster-steps 5 --min-seq-id 0.3
    """
}

clustDB.into {clustDB1; clustDB2; clustDB3 }


// Do the second clustering pass.
// Essentially, you align the cluster profiles against the profile consensus sequences
// then merge the clusters.
// This gets clusters down to ~10-20 % identity.
process mmseqsProfileClust {
    container "soedinglab/mmseqs2"

    input:
    set "clusters", "clusters.index" from clustDB2
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB3

    output:
    set "profile_clusters", "profile_clusters.index" into profileClustDB

    """
    mmseqs result2profile sequence sequence clusters profile

    mkdir -p tmp
    mmseqs search profile profile_consensus result tmp -s 7.5 -e 0.05 --add-self-matches --num-iterations 3
    mmseqs clust profile result profile_clusters_consensus

    mmseqs mergeclusters sequence profile_clusters clusters profile_clusters_consensus
    """
}


clustDB1.concat( profileClustDB ).into {allClusters1; allClusters2 }


// Extract some information and statistics about the clusters
process mmseqsExtractClusters {
    container "soedinglab/mmseqs2"

    publishDir "clusters"

    input:
    set file(clusters), file(clusters_index) from allClusters1
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB2

    output:
    file "${clusters}.tsv" into clustersTSV
    file "${clusters}_rep.fasta" into clustersRepFasta
    file "${clusters}_stats.tsv" into clustersStats

    """
    mmseqs createtsv \
      sequence \
      sequence \
      "${clusters}" \
      "${clusters}.tsv"

    mmseqs result2repseq \
      sequence \
      "${clusters}" \
      "${clusters}_rep"

    mmseqs result2flat \
      sequence \
      sequence \
      "${clusters}_rep" \
      "${clusters}_rep.fasta" \
      --use-fasta-header

    mmseqs align sequence sequence "${clusters}" align -a
    mmseqs convertalis sequence sequence align "${clusters}_stats.tsv"
    """
}


// Construct multiple sequence alignments from the clusters.
process getClusterMSAs {
    container "soedinglab/mmseqs2"
    publishDir "cluster_msas"

    input:
    set file(clusters), file(clusters_index) from allClusters2
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB3

    output:
    set "${clusters}_msa", "${clusters}_msa.index" into clustersMSAs
    file "${clusters}_msa.fasta" into clustersMSAFasta

    """
    mmseqs result2msa sequence sequence "${clusters}" "${clusters}_msa"
    mmseqs result2flat sequence sequence "${clusters}_msa" "${clusters}_msa.fasta"
    """

    /*
      mmseqs cluster DB DB_clu tmp
      mmseqs createseqfiledb DB DB_clu DB_clu_seq
      mmseqs apply DB_clu_seq DB_clu_seq_msa -- clustalo -i -  --threads=1
    */
}


process splitMSAs {
    publishDir "cluster_msas"

    input:
    file msas from clustersMSAFasta.filter { clu, idx -> clu.index}

    output:
    file "*.fasta" into indivFastas

    """
#!/usr/bin/env python3

current_name = None
current = []
last_was = False
with open("${msas}", "r") as handle:
    for line in handle:
        if line.startswith(">"):
            if last_was:
                if len(current) > 1:
                    with open(current_name, "w") as out_handle:
                        print("".join(current[:-1]), file=out_handle)
                current = []
                current_name = line.lstrip(">").strip().split(" ", 1)[0] + ".fasta"
            last_was = True
        else:
            last_was = False

        current.append(line)
    """
}


process estimateTrees {
    container "quay.io/biocontainers/fasttree:2.1.10--h470a237_2"
    publishDir "trees"

    input:
    file msa from indivFastas.flatMap()

    output:
    file "${msa.baseName}.nwk" into indivTrees

    """
    fasttree -fastest -quiet < "${msa}" > "${msa.baseName}.nwk"
    """
}

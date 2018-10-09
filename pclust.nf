#!/usr/bin/env nextflow

params.genomes = "$baseDir/data/*.{fasta,gff3}"

genomes = Channel.fromFilePairs( params.genomes, flat: true )

genomes.tap { genomes1 }

process extractProteins {
    container "quay.io/biocontainers/genometools-genometools:1.5.10--h470a237_1"

    input:
    set val(label), file(fasta), file(gff) from genomes1

    output:
    file "${label}.faa" into proteins

    """
    gt extractfeat \
      -type CDS \
      -translate \
      -matchdescstart \
      -join \
      -retainids \
      -seqfile "${fasta}" \
      "${gff}" \
      > "${label}.faa"
    """
}

proteins.into { proteins1; proteins2 }

// Combine the fasta files in preparation for clustering
// Need step to rename fastas to include filename.
process combineFasta {

    input:
    file "*.fa*" from proteins1.collect()

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
    publishDir "clusters"

    input:
    file fasta from combinedFasta

    output:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  into sequenceDB

    """
    mmseqs createdb "${fasta}" sequence --max-seq-len 14000
    """
}


sequenceDB.into {
    sequenceDB1; 
    sequenceDB2;
    sequenceDB3;
    sequenceDB4;
    sequenceDB5;
    sequenceDB6
    }

// Do the first clustering pass.
// This tends to get to about 50% identity.


process mmseqsClust {
    container "soedinglab/mmseqs2"
    publishDir "clusters"

    input:
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB1

    output:
    set "clusters", "clusters.index" into clustDB

    """
    mkdir -p tmp
    mmseqs cluster sequence clusters tmp -s 5 --cluster-steps 5 --min-seq-id 0.3
    """
}

clustDB.into {
    clustDB1;
    clustDB2;
    clustDB3;
    clustDB4;
    clustDB5
}

// Do the second clustering pass.
// Essentially, you align the cluster profiles against the profile consensus sequences
// then merge the clusters.
// This gets clusters down to ~10-20 % identity.

process mmseqsProfileClust {
    container "soedinglab/mmseqs2"
    publishDir "clusters"

    input:
    set "clusters", "clusters.index" from clustDB1
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup", "sequence_h", "sequence_h.index"  from sequenceDB2

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

profileClustDB.into { profileClustDB1; profileClustDB2 }


allClusters = Channel.value("identity").combine(clustDB2).concat(
    Channel.value("profile").combine( profileClustDB1 )
    )

allClusters.into { allClusters1; allClusters2 }

// Extract some information and statistics about the clusters
process mmseqsExtractClusters {
    container "soedinglab/mmseqs2"

    publishDir { "clusters/${type}" }
    tag { type }

    input:
    set val(type), file(clusters), file(clusters_index) from allClusters1
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
        "sequence_h", "sequence_h.index" from sequenceDB3

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
    publishDir "msas"
    tag { type }

    input:
    set val(type), file(clusters), file(clusters_index) from allClusters2
    set "sequence", "sequence.dbtype", "sequence.index", "sequence.lookup",
        "sequence_h", "sequence_h.index" from sequenceDB4

    output:
    set val(type), "${clusters}_msa", "${clusters}_msa.index" into clustersMSAs
    set val(type), "${clusters}_msa.fasta" into clustersMSAFasta

    """
    mmseqs result2msa sequence sequence "${clusters}" "${clusters}_msa"
    mmseqs result2flat sequence sequence "${clusters}_msa" "${clusters}_msa.fasta"
    """

    //  mmseqs cluster DB DB_clu tmp
    //  mmseqs createseqfiledb DB DB_clu DB_clu_seq
    //  mmseqs apply DB_clu_seq DB_clu_seq_msa -- clustalo -i -  --threads=1
}

/*

clustersMSAFasta.tap { clustersMSAFasta1 }

process splitMSAs {
    publishDir { "msas/${type}"}
    tag { type }

    input:
    set val(type), file(msas) from clustersMSAFasta1

    output:
    set val(type), file("*.fasta") into clustersMSAFastaIndiv

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
    set val(type), file(msa) from clustersMSAFastaIndiv.transpose()

    output:
    file "${msa.baseName}.nwk" into indivTrees

    """
    fasttree -fastest -quiet < "${msa}" > "${msa.baseName}.nwk"
    """
}

combinedFasta.tap { combinedFasta1 }

process effectorp {
    container "effectorp"
    publishDir "effectorp"

    input:
    file fasta from combinedFasta1

    output:
    file "table.tsv" into effectorpResults

    """
    EffectorP.py -s -i "${fasta}" > table.tsv
    """
}

*/

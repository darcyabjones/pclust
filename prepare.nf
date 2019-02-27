#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    =================================
    pclust/prepare
    =================================

    Prepares predicted proteomes for processing through the pclust pipeline(s).
    Because we use multiple proteomes in our studies, we need a way of making
    individual protein ids unique within a proteome and an easy way of finding
    which species a protein comes from.

    This pipeline adds tags to the start of sequences to make this possible.

    Usage:

    ```bash
    nextflow run -resume prepare.nf \
      --genomes "*.fasta" \
      --gffs "*.gff3" \
      --proteins "*.faa"
    ```

    ## Mandatory Arguments

    ```
    param                | description
    --------------------------------------------------
    `--genomes <fasta>`  | The input genome fasta files.

    `--gffs <gff3>`      | The gff3 files containing genes to extract
                         | from genomes.

    `--proteins <fasta>` | The protein fasta files to use directly.
    ```

    Note that either the genomes and gffs must be provided or the proteins
    must be provided. The pipeline can't run without input, or with only genomes
    and proteins (no gff).

    Outputs:

    * `sequences/gff_id_map.tsv`:
        The mapping of the original GFF3 ID protein names to the new ones.
    * `sequences/proteins.faa`:
        The combined processed fasta protein files.

    """.stripIndent()
}


params.genomes = false
params.gffs = false
params.proteins = false


if (params.help){
    helpMessage()
    exit 0
}


/*
 * Validate parameters.
 */
if ( params.gffs && !params.genomes ) {
    log.info "Extracting GFFs requires fasta genomes."
    exit 1
}

if ( !((params.gffs && params.genomes) || params.proteins) ) {
    log.info "You must input either gffs+genomes or proteins."
    exit 1
}


if ( params.genomes ) {
    genomes = Channel
        .fromPath( params.genomes )
        .map { [it.baseName, it] }
}


if ( params.gffs && params.genomes ) {

    gffs = Channel
        .fromPath( params.gffs )
        .map { [it.baseName, it] }

    /*
     * Duplicate the ID and Parent fields in the GFF into new attributes,
     * so that we can keep track of them after gt tidies them.
     */
    process duplicateIds {
        label "R"
        tag { label }

        input:
        set val(label), file("input.gff3") from gffs

        output:
        set val(label), file("${label}.gff3") into gffsDuplicatedIds

        """
        duplicate_gff_id_attribute.R input.gff3 > "${label}.gff3"
        """
    }

    /*
     * Tidy GFF3s so that genometools doesn't panic.
     * Note that we rename input files to avoid name clashes.
     */
    process tidyGFFs1 {
        label "genometools"
        tag { label }

        input:
        set val(label), file("input.gff3") from gffsDuplicatedIds

        output:
        set val(label), file("${label}.gff3") into tidiedGffs1

        """
        gt gff3 \
          -tidy \
          -sort \
          input.gff3 \
        > "${label}.gff3"
        """
    }

    /*
     * Add the filename to the IDs in the GFF, so that we know which
     * genome a protein came from.
     */
    process addFnameToId {
        label "R"
        tag { label }

        input:
        set val(label), file("${label}.gff3") from tidiedGffs1

        output:
        set val(label), file("out.gff3") into gffsWithFilenames
        set val(label), file("out.tsv") into gffFilenameMaps

        """
        add_filename_to_gff.R "${label}.gff3" out.gff3 out.tsv
        """
    }

    /*
     * Combine all map tables into a single file.
     */
    combinedGffMaps = gffFilenameMaps.collectFile(
        name: "gff_id_map.tsv",
        storeDir: "sequences",
        sort: "hash",
        newLine: true
    )

    /*
     * Tidy GFF3s so that genometools doesn't panic.
     * For whatever reason, a second tidy is needed to iron out the kinks for
     * gt.
     * Note that we rename input files to avoid name clashes.
     */
    process tidyGFFs2 {
        label "genometools"
        tag { label }

        input:
        set val(label), file("input.gff3") from gffsWithFilenames

        output:
        set val(label), file("${label}.gff3") into tidiedGffs2

        """
        gt gff3 \
          -tidy \
          -sort \
          -retainids \
          input.gff3 \
        > "${label}.gff3"
        """
    }

    /*
     * Join gff and genome channels together.
     */
    genomesGffs = genomes.join(tidiedGffs2, by: 0)

    /*
     * Extract protein sequences from genomes and GFF.
     */
    process extractProteins {
        label "genometools"
        tag { label }

        input:
        set val(label), file(fasta), file(gff) from genomesGffs

        output:
        set val(label), file("${label}.faa") into gffProteins

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

} else {
    gffProteins = Channel.empty()
}


if ( params.proteins ) {
    /*
     * Using the protein option we just need to add filenames to sequences.
     */

    raw_proteins = Channel
            .fromPath( params.proteins )
            .map { [it.baseName, it] }

    /*
     * Add filename to the sequence ids.
     */
    process addFnameToProtId {
        label "posix"
        tag { label }

        input:
        set val(label), file("input.fasta") from raw_proteins

        output:
        set val(label), file("${label}.faa") into proteinProteins

        """
        sed "s~>\\s*~>${label}_~g" < input.fasta > "${label}.faa"
        """
    }

} else {
    proteinProteins = Channel.empty()
}


proteins = gffProteins.concat( proteinProteins )

/*
 * Combine all proteins into a single fasta file.
 */
combinedFasta = proteins.map {l, f -> f} .collectFile(
    name: "proteins.faa",
    storeDir: "sequences",
    sort: "hash",
    newLine: true
)

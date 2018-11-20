#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=workq
#SBATCH --time=1-00:00:00
#SBATCH --account=y95
#SBATCH --mail-type=ALL
#SBATCH --mail-user=darcy.a.jones@postgrad.curtin.edu.au
#SBATCH --export=NONE

module load nextflow/18.10.1-bin
#      --seqs              description
#      --db
#      --enrich_seqs
#      --enrich_db
#
#    Options:
#      --trees
#      --nomsa
#      --nomsa_refine
#      --enrich_seqs
#      --enrich_db
#      --enrich_msa

nextflow run -resume -profile pawsey_zeus ./pclust.nf --max_cpus 28 --nomsa_refine --seqs "sequences/proteins.faa" --enrich_seqs "databases/uniclust50_2018_08/uniclust50_2018_08_consensus.fasta" --enrich_msa

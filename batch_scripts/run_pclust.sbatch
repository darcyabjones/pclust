#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=longq
#SBATCH --time=2-00:00:00
#SBATCH --account=y95
#SBATCH --mail-type=ALL
#SBATCH --mail-user=darcy.a.jones@postgrad.curtin.edu.au
#SBATCH --export=NONE

module load nextflow/19.01.0.5050-bin
module load singularity/3.3.0
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

nextflow run -resume \
  -profile pawsey_zeus,singularity \
  ./main.nf \
  --seqs nr_leq5000.fasta \
  --msas clusters/msa \
  --enrich_seqs enrich_seqs.fasta \
  --noremote \
  --outdir run

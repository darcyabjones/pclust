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

nextflow run -resume -profile pawsey_zeus ./phmm.nf --max_cpus 28 --nopfam --nouniref --msas "msas/mmseqs/*.fasta"

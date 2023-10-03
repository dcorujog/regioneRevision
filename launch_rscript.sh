#!/bin/bash
#SBATCH -J regioner
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem-per-cpu=3G
#SBATCH -o regioner_%j.out
#SBATCH -e regioner_%j.err
#SBATCH --time=3-00:00:00

module load R/4.2.2-foss-2021b

Rscript ${SLURM_SUBMIT_DIR}/permtest_hpc_genomic_features.R

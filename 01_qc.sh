#!/bin/bash
#SBATCH --job-name=fastqc_sample1
#SBATCH --output=fastqc_sample1.out
#SBATCH --error=fastqc_sample1.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --partition=pibu_el8

module load FastQC/0.11.9-Java-11
fastqc --threads $SLURM_CPUS_PER_TASK reads.sample1.fastq.gz
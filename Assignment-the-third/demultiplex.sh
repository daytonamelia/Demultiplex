#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=1-0
#SBATCH --job-name=demultiplex
#SBATCH --output=slurm_out/testing/slurm%j_blastp.out
#SBATCH --error=slurm_out/testing/slurm%j_blastp.err

read1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
read2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
read3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
read4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz

./demultiplex.py \
-r1 $read1 \
-r2 $read2 \
-r3 $read3 \
-r4 $read4 \
-i /projects/bgmp/shared/2017_sequencing/indexes.txt
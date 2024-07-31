#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=1-0
#SBATCH --job-name=demultiplex
#SBATCH --output=slurm_out/testing/slurm%j_blastp.out
#SBATCH --error=slurm_out/testing/slurm%j_blastp.err

read1=/projects/bgmp/adayton/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/testinput_R1_01.fastq.gz
read2=/projects/bgmp/adayton/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/testinput_R2_01.fastq.gz
read3=/projects/bgmp/adayton/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/testinput_R3_01.fastq.gz
read4=/projects/bgmp/adayton/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/testinput_R4_01.fastq.gz

./demultiplex.py \
-r1 $read1 \
-r2 $read2 \
-r3 $read3 \
-r4 $read4 \
-i /projects/bgmp/shared/2017_sequencing/indexes.txt
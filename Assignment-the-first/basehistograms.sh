#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=1-0
#SBATCH --job-name=demultiplex_histograms
#SBATCH --output=slurm_out/slurm%j_blastp.out
#SBATCH --error=slurm_out/slurm%j_blastp.err
#SBATCH --mail-user=daytonamelia@gmail.com
#SBATCH --mail-type=ALL

read1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
index1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
index2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
read2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz

/usr/bin/time -v ./perbasehist.py \
-f $read1 -l 101 -u True -y 42 \
-o ./R1_perbasehist.png

/usr/bin/time -v ./perbasehist.py \
-f $read2 -l 101 -u True -y 42 \
-o ./R2_perbasehist.png

/usr/bin/time -v ./perbasehist.py \
-f $index1 -l 8 -u True -y 42 \
-o ./index1_perbasehist.png

/usr/bin/time -v ./perbasehist.py \
-f $index2 -l 8 -u True -y 42 \
-o ./index2_perbasehist.png
# 2024 Bi622 Demultiplexing Lab Notebook

 ------------------
| Thursday 240725  |
 ------------------

## Assignment the First

From Assignment the first README.md: "Our goal is to look through a lane of sequencing generated from the 2017 BGMP cohort’s library preps and determine the level of index swapping and undetermined index-pairs, before and after quality filtering of index reads. In order to do this, we must first demultiplex the data. In Assignment the First, we will develop a strategy to de-multiplex samples to create 48 FASTQ files that contain acceptable index pairs (read1 and read2 for 24 different index pairs), two FASTQ files with index-hopped reads-pairs, and two FASTQ files undetermined (non-matching or low quality) index-pairs."

Starting Assignment the First!

## AtF: Part 1 - Quality Score Distribution per-nucleotide

Checking for files and sizes:

```
(base) [10:45:02 /projects/bgmp/shared/2017_sequencing] adayton:n0349.talapas.uoregon.edu
$ ls -lh
total 46G
-rw-r-xr--+ 1 coonrod  is.racs.pirg.bgmp  20G Jul 30  2018 1294_S1_L008_R1_001.fastq.gz*
-rw-r-xr--+ 1 coonrod  is.racs.pirg.bgmp 2.6G Jul 30  2018 1294_S1_L008_R2_001.fastq.gz*
-rw-r-xr--+ 1 coonrod  is.racs.pirg.bgmp 2.8G Jul 30  2018 1294_S1_L008_R3_001.fastq.gz*
-rw-r-xr--+ 1 coonrod  is.racs.pirg.bgmp  21G Jul 30  2018 1294_S1_L008_R4_001.fastq.gz*
drwxrws---+ 2 coonrod  is.racs.pirg.bgmp 8.0K Jul  1  2022 demultiplexed/
-rwxrwxr-x+ 1 sdwagner is.racs.pirg.bgmp  631 Aug  9  2021 indexes.txt*
-rw-r-xr--+ 1 coonrod  is.racs.pirg.bgmp  327 Aug 16  2017 README.txt*
```

Checking for line count in file 1:

```
(base) [10:29:06 /projects/bgmp/shared/2017_sequencing] adayton:n0349.talapas.uoregon.edu
$ zcat 1294_S1_L008_R1_001.fastq.gz | wc -l 
1452986940
```

So there are 363,246,735 records.

Head of each four files:

```
(base) [10:28:30 /projects/bgmp/shared/2017_sequencing] adayton:n0349.talapas.uoregon.edu
$ zcat 1294_S1_L008_R1_001.fastq.gz | head
@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1
GNCTGGCATTCCCAGAGACATCAGTACCCAGTTGGTTCAGACAGTTCCTCTATTGGTTGACAAGGTCTTCATTTCTAGTGATATCAACACGGTGTCTACAA
+
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ
@K00337:83:HJKJNBBXX:8:1101:1286:1191 1:N:0:1
CNACCTGTCCCCAGCTCACAGGACAGCACACCAAAGGCGGCAACCCACACCCAGTTTTACAGCCACACAGTGCCTTGTTTTACTTGAGGACCCCCCACTCC
+
A#AAFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJAJJJJJJJJJJJJJJFJJJJJFFFFJJJJJJJJJJJJJJJJJJ77F
@K00337:83:HJKJNBBXX:8:1101:1347:1191 1:N:0:1
GNGGTCTTCTACCTTTCTCTTCTTTTTTGGAGGAGTAGAATGTTGAGAGTCAGCAGTAGCCTCATCATCACTAGATGGCATTTCTTCTGAGCAAAACAGGT
(base) [10:37:44 /projects/bgmp/shared/2017_sequencing] adayton:n0349.talapas.uoregon.edu
$ zcat 1294_S1_L008_R2_001.fastq.gz | head
@K00337:83:HJKJNBBXX:8:1101:1265:1191 2:N:0:1
NCTTCGAC
+
#AA<FJJJ
@K00337:83:HJKJNBBXX:8:1101:1286:1191 2:N:0:1
NACAGCGA
+
#AAAFJJJ
@K00337:83:HJKJNBBXX:8:1101:1347:1191 2:N:0:1
NTCCTAAG
(base) [10:41:14 /projects/bgmp/shared/2017_sequencing] adayton:n0349.talapas.uoregon.edu
$ zcat 1294_S1_L008_R3_001.fastq.gz | head
@K00337:83:HJKJNBBXX:8:1101:1265:1191 3:N:0:1
NTCGAAGA
+
#AAAAJJF
@K00337:83:HJKJNBBXX:8:1101:1286:1191 3:N:0:1
NCGCTGTT
+
#AAAFJ-A
@K00337:83:HJKJNBBXX:8:1101:1347:1191 3:N:0:1
NTTAGGAC
(base) [10:41:26 /projects/bgmp/shared/2017_sequencing] adayton:n0349.talapas.uoregon.edu
$ zcat 1294_S1_L008_R4_001.fastq.gz | head
@K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1
NTTTTGATTTACCTTTCAGCCAATGAGAAGGCCGTTCATGCAGACTTTTTTAATGATTTTGAAGACCTTTTTGATGATGATGATGTCCAGTGAGGCCTCCC
+
#AAFAFJJ-----F---7-<FA-F<AFFA-JJJ77<FJFJFJJJJJJJJJJAFJFFAJJJJJJJJFJF7-AFFJJ7F7JFJJFJ7FFF--A<A7<-A-7--
@K00337:83:HJKJNBBXX:8:1101:1286:1191 4:N:0:1
NTGTGTAGACAAAAGTTTTCATGAGTCTGTAAGCTGTCTATTGTCTCCTGAAAAGAAACCAGAAGTTTTCCCCTAAATGTGTTTAGAATGCTTATTCTAAT
+
#A-AFFJJFJJJJJJJJJJJJJJJJ<JAJFJJJJF<JFJJJAJJJJJJJJJJJJJJJJJJJFJJJAJJFJJJFJJJF<JJA-JJJ-<AFAF--FF<JAFJF
@K00337:83:HJKJNBBXX:8:1101:1347:1191 4:N:0:1
NAAATGCCATCTAGTGATGATGAGGCTACTGCTGACTCTCAACATTCTACTCCTCCAAAAAAGAAGAGAAAGATTCCAACCCCCAGAACCGATGACCGGCA
(base) [10:41:31 /projects/bgmp/shared/2017_sequencing] adayton:n0349.talapas.uoregon.edu
$ 

```

Looks like R1 and R4 are reads, R2 and R3 are indexes.

R3 is the reverse complement of R2 - so these indexes match up. If they didnt match up, there was a hop. If one was not an index we recognize, then its unknown

I'm confused about the connections between the barcodes to R1 and R4.

OK I think R2 is the barcode for R1, R3 is the barcode for R4.

Looks like the R2 barcode is NOT reverse complemented to match our barcode list so I'm naming that index1.

SO
read 1: read1
read 2: index1
read 3: index2
read 4: read2

The reads appear to be 101 characters long just from copy and pasting some lines into a character counter. The indexes are 8 characters long.

```
(base) [11:33:36 /projects/bgmp/shared/2017_sequencing] adayton:n0349.talapas.uoregon.edu
$ zcat 1294_S1_L008_R1_001.fastq.gz | head -n 10000 | grep -c "a"
0
(base) [11:33:43 /projects/bgmp/shared/2017_sequencing] adayton:n0349.talapas.uoregon.edu
$ zcat 1294_S1_L008_R1_001.fastq.gz | head -n 10000 | grep -c "7"
3053
```

So it's +33 encoding as a "7" would only show up in +33 and an "a" would presumably show up in +64.

Are there Ns in every single index? Leslie said that if there are Ns we can disregard it. When looking at the amount of Ns I used:

```
$ zcat 1294_S1_L008_R2_001.fastq.gz | grep -A 1 "^@" | grep -v "^@" | grep -v "^--$" | grep -c "N"
3976613
```

There are 363246735 records, so 1.1% of the indexes in R2 have Ns.

```
$ zcat 1294_S1_L008_R3_001.fastq.gz | grep -A 1 "^@" | grep -v "^@" | grep -v "^--$" | grep -c "N"
3328051
```

0.92% of the indexes in R3 have Ns.
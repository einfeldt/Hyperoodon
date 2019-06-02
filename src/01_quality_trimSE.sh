#!/bin/sh

###########################################################
## 1 - Quality trim raw reads
###########################################################

# Detect number of CPUs available
#cpus="$(grep -c ^processor /proc/cpuinfo)"

# Trim reads
cd ./adapters
mkdir ../output/data_trimmed
for file in ../data_fastqgz/*R1_001.fastq.gz; do 
	filename="$(basename $file)"
	individualF="${filename%%.*}"
	# For Linux: TrimmomaticSE -threads $cpus -phred33 ../data_fastqgz/$individualF".fastq.gz" ../data_trimmed/$individualF".fastq.gz" LEADING:25 TRAILING:25 SLIDINGWINDOW:4:23 MINLEN:30 ILLUMINACLIP:TruSeq3-PE.fa:2:40:15
	java -jar /usr/local/bin/trimmomatic-0.33.jar SE -threads 4 -phred33 ../data_fastqgz/$individualF".fastq.gz" ../output/data_trimmed/$individualF".fastq.gz" LEADING:25 TRAILING:25 SLIDINGWINDOW:4:23 MINLEN:30 ILLUMINACLIP:TruSeq3-PE_all.fa:2:40:15
done

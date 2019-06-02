#!/bin/sh

###########################################################
## 2 - Filter mtDNA reads using reference
###########################################################

# Build bowtie index
cd ./ref_fasta
mkdir ../output/data_aligned
~/bowtie2-2.3.4.1-macos-x86_64/bowtie2-build ./*.fasta ./bt2_index --threads 4 --quiet

# Map reads to reference
for file in ../output/data_trimmed/*R1_001.fastq.gz; do 
	filename="$(basename $file)"
	individualF="${filename%%.*}"
	individualR="${individualF/R1_001/R2_001}"
	individual="${individualF%%_*}"
	# Align with Bowtie2, keep all reads aligning concordantly or nonconcordantly (non-concordant reads likely due to circular genome)
	~/bowtie2-2.3.4.1-macos-x86_64/bowtie2 -p 4 --very-sensitive-local -x ./bt2_index -q -1 ../output/data_trimmed/$individualF".fastq.gz" -2 ../output/data_trimmed/$individualR".fastq.gz" | samtools view -bSF4 - > ../output/data_aligned/$individual".bam"
	# Write to fastq
	samtools fastq ../output/data_aligned/$individual".bam" > ../output/data_aligned/$individual"_mappedreads.fastq"
done

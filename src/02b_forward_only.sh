#!/bin/sh

###########################################################
## 2b - Assemble mt genome using forward reads only
###########################################################

###########################################################
## 2b-2 - Filter mtDNA reads using reference
###########################################################

cd ./ref_fasta
mkdir ../output/data_aligned_b
~/bowtie2-2.3.4.1-macos-x86_64/bowtie2-build ./*.fasta ./bt2_index --threads 4 --quiet

# Map reads to reference
for file in ../output/data_trimmed_b/*R1_001.fastq.gz; do 
	filename="$(basename $file)"
	individualF="${filename%%.*}"
	individual="${individualF%%_*}"
	# Align with Bowtie2, keep all reads aligning concordantly or nonconcordantly (non-concordant reads likely due to circular genome)
	~/bowtie2-2.3.4.1-macos-x86_64/bowtie2 -p 4 --very-sensitive-local -x ./bt2_index -q -U ../output/data_trimmed_b/$individualF".fastq.gz" | samtools view -bSF4 - > ../output/data_aligned_b/$individual".bam"
	# Write to fastq
	samtools fastq ../output/data_aligned/$individual".bam" > ../output/data_aligned/$individual"_mappedreads.fastq"
done

## Unfinished, as neither of thee F-only reads files were high enough quality to map to mt genome



## TESTING

## Try MIRA
cd ../
mkdir ./output/data_finished
cd ./output/data_finished

ProjectName="$individual"

echo '# Paired Illumina data
project = '$individual'
job = genome,mapping,accurate

parameters = -NW:mrnl=0 -NW:cac=no -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no

readgroup
is_reference
data = '"../../ref_fasta/NBW_whole_mitochondrion_NC_005273_1.fasta"'
strain = REFERENCE

readgroup = reads
data = '"../data_trimmed_b/Hyam-6_S386_L001_R1_001.fastq.gz"'
technology = solexa
strain = SPECIES' > $individual"_guided_manifest.txt"

~/mira_4.9.6_darwin15.4.0_x86_64_static/bin/mira ./$individual"_guided_manifest.txt" >&$individual"_guided_log_assembly.txt"





## TESTING
# Convert fastq to fasta
fastq_to_fasta -i ../output/data_trimmed_b/$individualF".fastq" -o ../output/data_trimmed_b/$individualF".fasta" -Q33
blastn -subject ./NBW_whole_mitochondrion_NC_005273_1.fasta -query ../output/data_trimmed_b/$individualF".fasta" -evalue 2 > ../output/data_trimmed_b/$individual"_blastn_results_relaxed.txt"
cat ../output/data_trimmed_b/$individual"_blastn_results.txt" | grep " hits " | awk '{print $2}' | sort | uniq -c
cat ../output/data_trimmed_b/$individual"_blastn_results_relaxed.txt" | grep " hits " | awk '{print $2}' | sort | uniq -c

# Align to Megaptera mtDNA genome
cd ../ref_fasta_megaptera
~/bowtie2-2.3.4.1-macos-x86_64/bowtie2-build ./*.fasta ./bt2_index --threads 4 --quiet
~/bowtie2-2.3.4.1-macos-x86_64/bowtie2 -p 4 --very-sensitive-local -x ./bt2_index -q -U ../output/data_trimmed_b/$individualF".fastq.gz" | samtools view -bSF4 - > ../output/data_aligned_b/$individual".bam"
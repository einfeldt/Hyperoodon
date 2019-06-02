#!/bin/sh

###########################################################
## 3 - Assemble mtDNA genome
###########################################################

# Detect number of CPUs available
#cpus="$(grep -c ^processor /proc/cpuinfo)"

mkdir ./output/data_assembled
cd ./output/data_assembled

# Define reference
reffile=`ls ../../ref_fasta/*.fasta`
refname="$(basename $reffile)"

for file in ../data_aligned/*mappedreads.fastq; do 
	# Define variables
	filename="$(basename $file)"
	individual="${filename%_*}"	
	ProjectName="$individual"
	FReads=../data_aligned/$individual"-1.fastq"

# Write manifest file
echo '# Single-end Illumina data
project = '$individual'
job = genome,mapping,accurate

parameters = -NW:mrnl=0 -NW:cac=no -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no

readgroup
is_reference
data = '../../ref_fasta/$refname'
strain = REFERENCE

readgroup = reads
data = '../data_aligned/$file'
technology = solexa
strain = SPECIES' > $individual"_manifest.txt"

# Run MIRA
~/mira_4.9.6_darwin15.4.0_x86_64_static/bin/mira ./$individual"_manifest.txt" >&$individual"_log_assembly.txt"
done

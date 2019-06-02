#!/bin/sh

###########################################################
## mtExtract
###########################################################

## Purpose: assemble mtDNA from Illumina low-coverage genomic sequence data
# 1 - Quality trim all reads with Trimmomatic
# 2 - Map reads to reference with Bowtie2, keep only reads that map
# 3 - Reference-guided assembly using MIRA
# 4 - Assembly statistics
# 5 - De novo assemble mtDNA genomes with adequate coverage, remap

# Requires: 
# Trimmomatic
# Bowtie2
# samtools (process alignment output)
# MIRA assembler
# fastx toolkit (reverse complement)
# EMBOSS (merger)
# Pilon (Polish assembly using short reads)

# Output: assembly in .fasta and .caf format

## Detect number of CPUs available
#cpus="$(grep -c ^processor /proc/cpuinfo)"

## Run pipeline
cd ~/Desktop/mtExtract
mkdir ./output/
echo "Started quality trim:" `date` > ./output/timing_log.txt
bash ./src/01_quality_trimSE.sh
echo "Started filtering:" `date` >> ./output/timing_log.txt
bash ./src/02_filter_readsSE.sh
echo "Started assembling:" `date` >> ./output/timing_log.txt
bash ./src/03_assemble_mtDNASE.sh
echo "Gathering assembly statistics:" `date` >> ./output/timing_log.txt
bash ./src/04_statistics.sh
echo "Finished preliminary assembling on" `date` >> ./output/timing_log.txt
echo "Started de novo assembly on" `date` >> ./output/timing_log.txt
bash ./src/05_assemble_finishedSE.sh
echo "Finished denovo assembly on" `date` >> ./output/timing_log.txt
bash ./src/06_denovo_statistics.sh
echo "Finished writing denovo assembly statistics" `date` >> ./output/timing_log.txt
echo "Started polishing guided assemblies that were not successfully assembled de novo on" `date` >> ./output/timing_log.txt
bash ./src/07_polish_guided_assembliesSE.sh
echo "Finished extracting mitochondrial genomes on" `date` >> ./output/timing_log.txt

## Organize
cd ~/Desktop/mtExtract
# Write all outputs to single fasta
cat ./output/final_assemblies/*_final.fasta >> ./output/final_assemblies/01_all_samples.fasta
# Copy reference to output folder
cp -r ./ref_fasta ./output/reference
# Clean up reference folder
rm ./ref_fasta/bt2_index.*
rm ./ref_fasta/*.fasta_100bp
# Rename output by reference fasta ID
reffile=`ls ./ref_fasta/*.fasta`
refname="$(basename $reffile)"
refroot="${refname%.*}"
mv ./output ./output_$refroot


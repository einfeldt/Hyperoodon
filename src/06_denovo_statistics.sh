#!/bin/sh

###########################################################
## 6 - Denovo assembly statistics
###########################################################

# Detect number of CPUs available
#cpus="$(grep -c ^processor /proc/cpuinfo)"

# Define reference
reffile=`ls ./ref_fasta/*.fasta`
referenceBases=$(tail -n +2 $reffile | wc -c)

## Create summary file
echo "sample,totalReads,trimmedReads,mappedReads,assembledReads,averageCoverage,maxCoverage,consensusQuality,assemblyLength,unmappedBases,mappedBases,finalAssemblyBases" > ./output/"denovo_assembly_stats.csv"

for file in `ls -d ./output/data_finished/*_remap_assembly/`; do
	# Define variables
	dirbase="$(basename $file)"
	individual="${dirbase%%_*}"
	seqfile="./data_fastqgz/"$individual"*R1_001.fastq.gz"
	filename="$(basename $seqfile)"	

	# Number of reads (R1 and R2 of paired reads)
	totalReads=$(gunzip -c ./data_fastqgz/$filename | wc -l | awk '{print $1/2}')
	# Number of trimmed reads (R1 and R2 of paired reads)
	trimmedReads=$(gunzip -c ./output/data_trimmed/$filename | wc -l | awk '{print $1/2}')
	# Number of mapped reads (=npaired*2)
	mappedReads=$(cat ./output/data_aligned/$individual"_mappedreads.fastq" | wc -l | awk '{print $1/4}')
	# Number of reads assembled (=npaired*2)
	assembledReads=$(expr `head -9 ./output/data_finished/$individual"_remap_assembly"/$individual"_remap_d_info/"$individual"_remap_info_assembly.txt" | tail -1 | awk '{print $NF}'` - 1)
	# Average coverage
	averageCoverage=$(head -15 ./output/data_finished/$individual"_remap_assembly"/$individual"_remap_d_info/"$individual"_remap_info_assembly.txt" | tail -1 | awk '{print $NF}')
	# Maximum coverage
	maxCoverage=$(head -95 ./output/data_finished/$individual"_remap_assembly"/$individual"_remap_d_info/"$individual"_remap_info_assembly.txt" | tail -1 | awk '{print $NF}')
	# Average consensus quality
	consensusQuality=$(head -100 ./output/data_finished/$individual"_remap_assembly"/$individual"_remap_d_info/"$individual"_remap_info_assembly.txt" | tail -1 | awk '{print $NF}')
	# Assembly length
	assemblyLength=$(head -79 ./output/data_finished/$individual"_remap_assembly"/$individual"_remap_d_info/"$individual"_remap_info_assembly.txt" | tail -1 | awk '{print $NF}')
	# Number basepairs in reference not covered by mapped reads
	unmappedSequence=$(tail -n +2 ./output/data_finished/$individual"_remap_assembly"/$individual"_remap_d_results/"$individual"_remap_out_SPECIES.unpadded.fasta" | tr -d '[:space:]')
	unmappedBases=$(grep -o "X" <<<"$unmappedSequence" | wc -l)
	# Number basepairs in reference covered by mapped reads
	totalBases=$(tail -n +2 ./output/data_finished/$individual"_remap_assembly"/$individual"_remap_d_results/"$individual"_remap_out_SPECIES.unpadded.fasta" | tr -d '[:space:]'| wc -c)
	mappedBases=$(expr $totalBases - $unmappedBases)
	# Number basepairs in reference covered by mapped reads
	finalAssemblyBases=$(tail -n +2 ./output/final_assemblies/$individual"_final.fasta" | tr -d '[:space:]'| wc -c)

	# Write variables to file
	echo $individual,$totalReads,$trimmedReads,$mappedReads,$assembledReads,$averageCoverage,$maxCoverage,$consensusQuality,$assemblyLength,$unmappedBases,$mappedBases,$finalAssemblyBases >> ./output/"denovo_assembly_stats.csv"
done

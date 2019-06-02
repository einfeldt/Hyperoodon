#!/bin/sh

###########################################################
## 7 - Polish guided mtDNA genomes for samples failing denovo assembly
###########################################################

## Define reference
reffile=`ls ./ref_fasta/*.fasta`
refname="$(basename $reffile)"

# Length of reference mtDNA genome (based on most common result from denovo assembly)
refLength=$(awk 'NF>1{print $NF}' ./output/denovo_assembly_stats.csv | sort | uniq -c | sort -nr | head -1 | awk '{print $2}')

# Conditionally overwrite files that do not match reference length with reference-guided assembly
for file in `ls -d ./output/data_finished/*_remap_assembly/`; do 
	filename="$(basename $file)"
	individual="${filename%%_*}"	
	finalAssemblyBases=$(tail -n +2 ./output/final_assemblies/$individual"_final.fasta" | tr -d '[:space:]'| wc -c)
	if [[ "$finalAssemblyBases" -eq "$refLength" ]]; then 	
		echo "Number of bases in assembly ("$finalAssemblyBases") matches reference length of "$refLength" bases. Ignoring file."
	else
		echo "Number of bases in assembly ("$finalAssemblyBases") does not match reference length of "$refLength" bases. Overwriting final assembly file with polished guided assembly."
		# Build reference from denovo mtDNA genome
		~/bowtie2-2.3.4.1-macos-x86_64/bowtie2-build ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES.unpadded.fasta" ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_bt2_index" --threads 4 --quiet
		# Align start of reference to de novo assembly with bowtie2
		~/bowtie2-2.3.4.1-macos-x86_64/bowtie2 -p 4 --very-sensitive-local -x ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_bt2_index" -f ./ref_fasta/$refname"_100bp" | samtools view -SF4 - > ./output/data_assembled/$individual"_assembly"/$individual".sam"
		# Conditionally reverse complement sequence to match reference orientation
		samFlag=$(awk '{print $2}' ./output/data_assembled/$individual"_assembly"/$individual".sam")
		if [[ $samFlag -eq 16 ]]; then 
			mv ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES.unpadded.fasta" ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES_original.unpadded.fasta"
			# Collapse to single line
			awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES_original.unpadded.fasta" | tail -n +2 > ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES_original_1line.unpadded.fasta"
			# Sanitize sequence by changing all non-ATGC characters to N
			cat ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES_original_1line.unpadded.fasta" | sed '/^[^>]/ s/[^AGTC]/N/g' > ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES_original_1line_sanitized.unpadded.fasta"
			# Reverse complement
			fastx_reverse_complement -i ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES_original_1line_sanitized.unpadded.fasta" -o ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES.unpadded.fasta"	
		fi
		
		## Polish assembly by remapping short reads
		~/bowtie2-2.3.4.1-macos-x86_64/bowtie2-build ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES.unpadded.fasta" ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_guided_bt2_index" --threads 4 --quiet
		# Align start of reference to de novo assembly with bowtie2
		~/bowtie2-2.3.4.1-macos-x86_64/bowtie2 -p 4 --very-sensitive-local -x ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_guided_bt2_index" -q --interleaved ./output/data_aligned/$individual"_mappedreads.fastq" | samtools view -bS | samtools sort - > ./output/data_assembled/$individual"_assembly"/$individual"_guided.bam"
		samtools index ./output/data_assembled/$individual"_assembly"/$individual"_guided.bam" ./output/data_assembled/$individual"_assembly"/$individual"_guided.bai"
		# Polish with Pilon (note: higher memory required if higher number of aligned sequences, use: -Xmx20G on more powerful machine)
		java -jar ~/pilon-1.22.jar --genome ./output/data_assembled/$individual"_assembly"/$individual"_d_results"/$individual"_out_SPECIES.unpadded.fasta" --frags ./output/data_assembled/$individual"_assembly"/$individual"_guided.bam" --output ./output/data_assembled/$individual"_assembly"/$individual"_final_tmp"
		# Rename fasta entry
		fasta_ID=">"$individual"_guided_assembled_pilon"
		cat ./output/data_assembled/$individual"_assembly"/$individual"_final_tmp.fasta" | sed "1s/.*/$fasta_ID/" > ./output/data_assembled/$individual"_assembly"/$individual"_final.fasta"
		
		## Re-map and assemble with MIRA to produce coverage information (using strict sequence length)
		# Write manifest file
echo '# Paired Illumina data
project = '$individual"_remap"'
job = genome,mapping,accurate

parameters = -NW:mrnl=0 -NW:cac=no -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no

readgroup
is_reference
data = './$individual"_assembly"/$individual"_final.fasta"'
strain = REFERENCE

readgroup = reads
data = '../data_aligned/$individual"_mappedreads.fastq"'
technology = solexa
strain = SPECIES' > ./output/data_assembled/$individual"_remap_manifest.txt"
		# Run MIRA
		cd ./output/data_assembled
		~/mira_4.9.6_darwin15.4.0_x86_64_static/bin/mira ./$individual"_remap_manifest.txt" >&./$individual"_remap_log_assembly.txt"
		cd ../../
		
		## Copy .fasta and .caf files to final assembly folder
		rm ./output/final_assemblies/$individual"_out.caf"
		cp ./output/data_assembled/$individual"_remap_assembly"/$individual"_remap_d_results"/$individual"_remap_out.caf" ./output/final_assemblies/$individual"_out_guided.caf"
		rm ./output/final_assemblies/$individual"_final.fasta"
		cp ./output/data_assembled/$individual"_assembly"/$individual"_final.fasta" ./output/final_assemblies/$individual"_guided_final.fasta"
	fi
done

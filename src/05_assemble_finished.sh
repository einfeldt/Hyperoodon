#!/bin/sh

###########################################################
## 5 - Assemble finished mtDNA genome
###########################################################

# Detect number of CPUs available
#cpus="$(grep -c ^processor /proc/cpuinfo)"

# For assemblies with 100% coverage, assembly de novo then detect end overlap and re-map

# Define minimum thresholds for de novo assembly
minCoverage="15"
maxUnmapped="3"

mkdir ./output/data_finished
mkdir ./output/final_assemblies
cd ./output/data_finished

## Define reference
reffile=`ls ../../ref_fasta/*.fasta`
refname="$(basename $reffile)"

######################################################################################################################
## Assemble mtDNA genome de novo if sufficient coverage and lack of gaps in reference-guided assembly
for file in ../data_aligned/*mappedreads.fastq; do 

	# Define variables
	filename="$(basename $file)"
	individual="${filename%_*}"	
	ProjectName="$individual"
	unmappedSequence=$(tail -n +2 ../data_assembled/$individual"_assembly"/$individual"_d_results/"$individual"_out_SPECIES.unpadded.fasta")
	unmappedBases=$(grep -o "X" <<<"$unmappedSequence" | wc -l)
	averageCoverage=$(head -15 ../data_assembled/$individual"_assembly"/$individual"_d_info/"$individual"_info_assembly.txt" | tail -1 | awk '{print $NF}')
	averageCoverageInteger=${averageCoverage%.*}	

	# Conditionally assemble de novo
	if (( "$unmappedBases" > "$maxUnmapped" )) || (( "$averageCoverageInteger" < "$minCoverage" )); then 
		echo "unmapped bases ="$unmappedBases", average coverage = "$averageCoverageInteger", did not assemble de novo: >"$unmappedBases" bases from reference unmapped or < "$minCoverage" reads average coverage in reference-guided assembly" > $individual"_denovo_log_assembly.txt"
	else
# Write manifest file
echo '# Paired Illumina data
project = '$individual'
job = genome,denovo,accurate

parameters = -NW:mrnl=0 -NW:cac=no -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no

readgroup = reads
data = '$file'
technology = solexa
strain = SPECIES' > $individual"_denovo_manifest.txt"

# Run MIRA
~/mira_4.9.6_darwin15.4.0_x86_64_static/bin/mira ./$individual"_denovo_manifest.txt" >&$individual"_denovo_log_assembly.txt"

#############################################################
## Split consensus based on start of reference provided and merge overlap

# Define reference sequences
refseq=$(tail -n +2 $reffile | tr -d '[:space:]')
refseq100=${refseq:0:100}
echo ">First 100bp of "$refname"
"$refseq100 > ../../ref_fasta/$refname"_100bp"
# Build reference from denovo mtDNA genome
~/bowtie2-2.3.4.1-macos-x86_64/bowtie2-build ./$individual"_assembly"/$individual"_d_results"/$individual"_LargeContigs_out_SPECIES.unpadded.fasta" ./$individual"_assembly"/$individual"_d_results"/$individual"_bt2_index" --threads 4 --quiet
# Align start of reference to de novo assembly with bowtie2
~/bowtie2-2.3.4.1-macos-x86_64/bowtie2 -p 4 --very-sensitive-local -x ./$individual"_assembly"/$individual"_d_results"/$individual"_bt2_index" -f ../../ref_fasta/$refname"_100bp" | samtools view -SF4 - > ./$individual"_assembly"/$individual".sam"
# Conditionally reverse complement sequence to match reference orientation
samFlag=$(awk '{print $2}' ./$individual"_assembly"/$individual".sam")
if [[ $samFlag -eq 16 ]]; then 
	mv ./$individual"_assembly"/$individual"_d_results"/$individual"_LargeContigs_out_SPECIES.unpadded.fasta" ./$individual"_assembly"/$individual"_d_results"/$individual"_LargeContigs_out_SPECIES_original.unpadded.fasta"
	# Collapse to single line
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ./$individual"_assembly"/$individual"_d_results"/$individual"_LargeContigs_out_SPECIES_original.unpadded.fasta" | tail -n +2 > ./$individual"_assembly"/$individual"_d_results"/$individual"_LargeContigs_out_SPECIES_original_1line.unpadded.fasta"
	# Sanitize sequence by changing all non-ATGC characters to N
	cat ./$individual"_assembly"/$individual"_d_results"/$individual"_LargeContigs_out_SPECIES_original_1line.unpadded.fasta" | sed '/^[^>]/ s/[^AGTC]/N/g' > ./$individual"_assembly"/$individual"_d_results"/$individual"_LargeContigs_out_SPECIES_original_1line_sanitized.unpadded.fasta"
	# Reverse complement
	fastx_reverse_complement -i ./$individual"_assembly"/$individual"_d_results"/$individual"_LargeContigs_out_SPECIES_original_1line_sanitized.unpadded.fasta" -o ./$individual"_assembly"/$individual"_d_results"/$individual"_LargeContigs_out_SPECIES.unpadded.fasta"	
fi
# Get sequence of largest contig assembled denovo
fastaseq=$(tail -n +2 ./$individual"_assembly"/$individual"_d_results/"$individual"_LargeContigs_out_SPECIES.unpadded.fasta" | tr -d '[:space:]' | cut -f1 -d">")
# Basepair in denovo assembly matching start of reference
firstbp_r2tmp=$(blastn -query ../../ref_fasta/$refname"_100bp" -subject ./$individual"_assembly"/$individual"_d_results"/$individual"_LargeContigs_out_SPECIES.unpadded.fasta" -outfmt "6 sstart")
firstbp_r2=$(expr $firstbp_r2tmp - 1)
lastbp_r1=$(expr $firstbp_r2)
lastbp_r2=$(echo $fastaseq | wc -m)
# Split sequences
fastaseq_1=${fastaseq:0:$lastbp_r1}
fastaseq_2=${fastaseq:$firstbp_r2:$lastbp_r2}
# Write sequences
echo ">"$individual"_front_sequence
"$fastaseq_2 > ./$individual"_assembly"/$individual"_seqfront.fasta" # end of this sequence will overlap with beginning of other
echo ">"$individual"_back_sequence
"$fastaseq_1 > ./$individual"_assembly"/$individual"_seqback.fasta" # beginning of this sequence will overlap with end of other
# Merge overlapping sequences with EMBOSS merger
~/EMBOSS-6.6.0/emboss/merger ./$individual"_assembly"/$individual"_seqfront.fasta" ./$individual"_assembly"/$individual"_seqback.fasta" -outfile ./$individual"_assembly"/$individual"_merge.txt" -outseq ./$individual"_assembly"/$individual"_merge.fasta"

# Linearize and change to caps
mkdir ./$individual"_final"
prelimSequence=$(tail -n +2 ./$individual"_assembly"/$individual"_merge.fasta" | tr -d '[:space:]')
finalSequence=$(echo $prelimSequence | awk '{print toupper($0)}')
echo ">"$individual"_denovo_assembled_splitmerged
"$finalSequence > ./$individual"_final"/$individual"_splitmerged.fasta"

###########################################################
## Polish assembly by remapping short reads

# Align with Bowtie2
# Build reference from denovo mtDNA genome
~/bowtie2-2.3.4.1-macos-x86_64/bowtie2-build ./$individual"_final"/$individual"_splitmerged.fasta" ./$individual"_final"/$individual"_splitmerged_bt2_index" --threads 4 --quiet
# Align start of reference to de novo assembly with bowtie2
~/bowtie2-2.3.4.1-macos-x86_64/bowtie2 -p 4 --very-sensitive-local -x ./$individual"_final"/$individual"_splitmerged_bt2_index" -q --interleaved ../data_aligned/$individual"_mappedreads.fastq" | samtools view -bS | samtools sort - > ./$individual"_final"/$individual"_splitmerged.bam"
samtools index ./$individual"_final"/$individual"_splitmerged.bam" ./$individual"_final"/$individual"_splitmerged.bai"
# Polish with Pilon
java -jar ~/pilon-1.22.jar --genome ./$individual"_final"/$individual"_splitmerged.fasta" --frags ./$individual"_final"/$individual"_splitmerged.bam" --output ./$individual"_final"/$individual"_final"

###########################################################
## Re-map and assemble with MIRA to produce coverage information (using strict sequence length)

# Write manifest file
echo '# Paired Illumina data
project = '$individual"_remap"'
job = genome,mapping,accurate

parameters = -NW:mrnl=0 -NW:cac=no -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no

readgroup
is_reference
data = './$individual"_final"/$individual"_final.fasta"'
strain = REFERENCE

readgroup = reads
data = '$file'
technology = solexa
strain = SPECIES' > $individual"_remap_manifest.txt"
# Run MIRA
~/mira_4.9.6_darwin15.4.0_x86_64_static/bin/mira ./$individual"_remap_manifest.txt" >&$individual"_remap_log_assembly.txt"

###########################################################
## Copy .fasta and .caf files to final assembly folder
cp ./$individual"_remap_assembly"/$individual"_remap_d_results"/$individual"_remap_out.caf" ../final_assemblies/$individual"_out.caf"
cp ./$individual"_final"/$individual"_final.fasta" ../final_assemblies/$individual"_final.fasta"

	fi
done
######################################################################################################################

# Runs mixcr on paired end RNA seq fastq files
# Takes 1 arguments: species "hsa, mmu, etc."
# Assumes files and names were generated using the weldr protocol
# Assumes fastqs were aligned to decontaminome and unaligned files were named decon_unal_R1(2).fastq.gz
# Ex: bash mixcr_PE_rnaseq_CDR3.sh hsa


if [[ ${1} == "" ]]; then 
echo "Please input either hsa or mmu as argument."
exit
fi

for sample in $(cat samples.txt); do


echo '#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=mixcr_PE_kaligner2_'${sample}'
#SBATCH --account=PAS1555
#SBATCH --output=mixcr_PE_kaligner2_'${sample}'.log

date


cd $SLURM_SUBMIT_DIR
if [ ! -d output/'${sample}' ]; then
	mkdir output/'${sample}'
fi
cd output/'${sample}'
pwd

# Use trimmed fastqs as data input
if [ -d assembly ]; then
	FQ1=$(ls assembly/S0_R1.t*.fastq.gz)
	FQ2=$(ls assembly/S0_R2.t*.fastq.gz)
	STORE="mixcr_umi_collapsed"
elif [ -f decon_unal_R1.fastq.gz ]; then
	FQ1=$(ls decon_unal_R1.fastq.gz)
	FQ2=$(ls decon_unal_R2.fastq.gz)
	STORE="mixcr"
elif [ -f trim_R1.fastq.gz ]; then
	FQ1=$(ls trim_R1.fastq.gz)
	FQ2=$(ls trim_R2.fastq.gz)
	STORE="mixcr"
else
	FQ1=$(ls ../../raw/'${sample}'_R1.fastq.gz)
	FQ2=$(ls ../../raw/'${sample}'_R2.fastq.gz)	
	STORE="mixcr"
fi

echo $FQ1
echo $FQ2

echo "MIXCR ALIGNMENT TO:" '${1}'
#mixcr align -f -p rna-seq -s '${1}' -OallowPartialAlignments=true -OsaveOriginalReads=true -OallowNoCDR3PartAlignments=true $FQ1 $FQ2 alignments.vdjca
mixcr align -f -p kAligner2 -s '${1}' -OallowPartialAlignments=true -OsaveOriginalReads=true -OallowNoCDR3PartAlignments=true --report report.txt $FQ1 $FQ2 alignments.vdjca

echo "ASSEMBLING PARTIAL SEQUENCES"
mixcr assemblePartial -f --report report.txt alignments.vdjca alignmentsRescued_1.vdjca

mixcr assemblePartial -f --report report.txt alignmentsRescued_1.vdjca alignmentsRescued_2.vdjca

echo "EXTENDING ALIGNMENTS"
mixcr extend --report report.txt -f alignmentsRescued_2.vdjca alignmentsRescued_2_extended.vdjca

echo "ASSEMBLING ALIGNMENTS"
#mixcr assemble -f alignmentsRescued_2_extended.vdjca clones.clns
mixcr assemble --report report.txt -f -a alignmentsRescued_2_extended.vdjca clones.clna

echo "EXPORTING CLONES"
#mixcr exportClones -f -o -t clones.clns clones.txt
mixcr exportClones -f -o -t clones.clna clones.txt
mixcr exportClones -f clones.clna clones_all.txt

echo "EXPORTING ALIGNMENTS"
#mixcr exportAlignmentsPretty alignments.vdjca > alignments.pretty
#mixcr exportAlignments -f --preset full -descrsR1 -descrsR2 alignments.vdjca alignments.txt
#mixcr exportAlignments -f --preset full -descrsR1 -descrsR2 alignmentsRescued_2_extended.vdjca alignments_extended.txt
mixcr exportAlignments -f --preset full -descrsR1 -descrsR2 -cloneIdWithMappingType clones.clna alignments_clones.txt

awk -F '\t' '"'"'{gsub(/\.0/, "", $2)}1'"'"' OFS='\t' clones.txt > tmp && mv tmp clones.txt

echo "Assembling contigs"
mixcr assembleContigs --report report.txt clones.clna full_clones.clns
mixcr exportClones -p fullImputed full_clones.clns full_clones.txt

mkdir ${STORE}_kaligner2
mv alignments_extended.txt alignments.txt alignments_clones.txt alignments.vdjca alignmentsRescued_1.vdjca \
alignmentsRescued_2.vdjca alignmentsRescued_2_extended.vdjca clones.clna clones.txt clones_all.txt full_clones.clns full_clones.txt report.txt \
${STORE}_kaligner2/

date
echo "Workflow complete"
' | sbatch

done

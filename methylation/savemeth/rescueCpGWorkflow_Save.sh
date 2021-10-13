### This script is a wordflow to resuce methylation from multimapped and unaligned reads
### Run this script in a directory containing all fastq files
### Run after rescueCpGWorkflow_Align.sh
	# There should exist ../bismark_SEalign and ../bismark_multimapLOC
	# These directories should contain SE alignment and multimapping location bam files
### Usage: bash rescueCpGWorkflow_Save.sh [sample list] [reference fasta genome]
	### Ex: bash rescueCpGWorkflow_Save.sh samples.txt hg38.fa
	### samples.txt contains sample name identifiers on each new line

# Outputs in a directory called savedMeth/

# Please see dependancies as listed in rescue_ambigMeth.py

set -e 

if [[ ${1} == "" ]]; then
	echo "Please supply a text file containing sample name identifiers on each new line"
	echo "Usage: bash rescueCpGWorkflow_Save.sh [sample list] [reference fasta genome]"
	exit
fi

if [[ ${2} == "" ]]; then
	echo "Please supply a fasta reference genome"
	echo "Usage: bash rescueCpGWorkflow_Save.sh [sample list] [reference fasta genome]"
	exit
fi

if [ ! -d ../bismark_SEalign ]; then
        echo "../bismark_SEalign doesn't exist ... Please run rescueCpGWorkflow_Align.sh first!"
	exit
fi

if [ ! -d ../bismark_multimapLOC ]; then
        echo "../bismark_multimapLOC doesn't exist ... Please run rescueCpGWorkflow_Align.sh first!"
	exit
fi

if [[ $(ls ../bismark_SEalign) == "" ]]; then
        echo "../bismark_SEalign is empty ... Please run rescueCpGWorkflow_Align.sh first!"
	exit	
fi

if [[ $(ls ../bismark_multimapLOC) == "" ]]; then
        echo "../bismark_multimapLOC is empty ... Please run rescueCpGWorkflow_Align.sh first!"
	exit	
fi

if [ ! -d savedMeth ]; then
        mkdir -p savedMeth
fi

for NAME in $(cat $1); do

	FASTQ=$(ls ${NAME}*_V1B_comb_SE_reads.NOTuniquemap.fq.gz)
	SEALIGN=$(ls ../bismark_SEalign/${NAME}*_V1B_comb_SE_reads.NOTuniquemap_bismark_bt2.bam)
	MULTIPOS=$(ls ../bismark_multimapLOC/${NAME}*_V1B_comb_SE_reads.NOTuniquemap.fq.gz_bismark_bt2.ambig.bam)
	REF=${2}
	OUTNAME=${SEALIGN#../bismark_SEalign/}
	OUTNAME=${OUTNAME%%_*}_FINALsaved.permeth
	AMBIGNAME=${NAME}_ambig_names.txt
	UNNAME=${NAME}_unaligned_names.txt
	echo '
	#PBS -l walltime=1:00:00
	#PBS -l nodes=1:ppn=28
	#PBS -l mem=112GB
	#PBS -N savemeth_rescue_'${NAME}'
	#PBS -S /bin/bash
	#PBS -j oe
	#PBS -W umask=002
	#PBS -A PAS1359


	date

	cd $PBS_O_WORKDIR
	cd savedMeth
	
	pwd

	source activate /users/PAS0472/osu9900/anaconda2/envs/myenv2.7

	echo "/users/PAS0472/osu9900/RRBS_scripts/savemeth/rescue_ambigMeth_v2.py -b" '../${SEALIGN}' "-f" '../${REF}' "-q" '../${FASTQ}' "-a" '../${MULTIPOS}' \
	"-o" '${OUTNAME}' "-anames" '../${AMBIGNAME}' "-unnames" '../${UNNAME}'

	echo

	python /users/PAS0472/osu9900/RRBS_scripts/savemeth/rescue_ambigMeth_v2.py -b ../'${SEALIGN}' -f ../'${REF}' -q ../'${FASTQ}' -a ../'${MULTIPOS}' \
	-o '${OUTNAME}' -anames '${AMBIGNAME}' -unnames '${UNNAME}'

	date
	echo "Workflow complete"
	' | qsub

done

# Forked 01/09/2020 by Altan Turkoglu (turkoglu.12@osu.edu) for use
# in Kirby Lab technical paper. Original author: Alex Pan (alexander.pan@osumc.edu)

# Demultiplex and align 10X samples
# Usage: bash 10X-workflow.sh [REFNAME] [MKFASTQ] [SAMPLESHEET.CSV] [OPTIONAL: VDJ.txt]
	# [REFNAME] MM or HSA
	# [MKFASTQ] cellrange mkfastq output directory OR directory containing fastqs
		# Fastqs MUST be named like so: SAMPLENAME_S1_L001_R1_001.fastq.gz
	# [SAMPLESHEET.CSV] A list of SAMPLENAMEs to run
	# [VDJ.txt] By default, this script runs all samples through cellranger counts. Please give a text file with each sample name on a new line that specifies the sample names to run cellranger vdj on. This option can be left blank.

REFNAME=${1}
MKFASTQ=${2}
SAMPLESHEET=${3}
ACCOUNT='PAS0472'

if [[ ${REFNAME} == "HSA" ]]; then
REF="/users/PAS0472/osu9900/anaconda2/bin/cellranger-3.1.0/refdata-cellranger-GRCh38-3.0.0"
VDJREF="/users/PAS0472/osu9900/anaconda2/bin/cellranger-3.1.0/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0"
elif [[ ${REFNAME} == "MM" ]]; then
REF="/users/PAS0472/osu9900/anaconda2/bin/cellranger-3.1.0/refdata-cellranger-mm10-3.0.0"
VDJREF=""
else
echo "Please input HSA or MM as first argument"
exit 1
fi

if [[ -d ${MKFASTQ%/}/outs/fastq_path ]]; then
MKFASTQ=$(ls -d ${MKFASTQ%/}/outs/fastq_path)
fi

if [[ ! -f ${SAMPLESHEET} ]]; then
echo "Please run this script in the directqory that contains ${SAMPLESHEET}"
exit
fi

if [[ ${4} == "" ]]; then 
echo "All samples will be run through cellranger count."
VDJFILE="NONE"
else
echo "The following will be run through cellranger vdj: "
cat ${4}
VDJFILE=${4}
fi

if [[ ${VDJFILE} == "NONE" ]]; then 
	COUNTLIST=$(awk -F ',' '{print $1}' ${SAMPLESHEET} | \sort | uniq)
else
	COUNTLIST=$(awk -F ',' '{print $1}' ${SAMPLESHEET} | grep -v -w -f ${VDJFILE} | \sort | uniq)
	VDJLIST=$(awk -F ',' '{print $1}' ${SAMPLESHEET} | grep -w -f ${VDJFILE} | \sort | uniq)
fi

for SAMPLE in $COUNTLIST; do

echo '#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name=count_'${SAMPLE}'
#SBATCH --account='${ACCOUNT}'

cd $PBS_O_WORKDIR
pwd

source activate ~osu9900/anaconda2/envs/myenv2.7/

echo '${SAMPLE}' '${REF}' '${MKFASTQ}'
cellranger count --id=counts_'${SAMPLE}' --transcriptome='${REF}'\
	--fastqs='${MKFASTQ}' --sample='${SAMPLE}'


' | sbatch

done

for SAMPLE in $VDJLIST; do

echo '#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name=vdj_'${SAMPLE}'
#SBATCH --account='${ACCOUNT}'

cd $PBS_O_WORKDIR
pwd
source activate ~osu9900/anaconda2/envs/myenv2.7/

echo '$SAMPLE' '${VDJREF}' '${MKFASTQ}'

cellranger vdj --id=vdj_'${SAMPLE}' --fastqs='${MKFASTQ}' --sample='${SAMPLE}' --reference='${VDJREF}'

' | sbatch

done



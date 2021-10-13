# Demultiplex and align 10X samples
# Usage: bash 10X-workflow.sh [BCL DIRECTORY] [SAMPLESHEET.CSV]
	# [BCL DIRECTORY] Unpacked tar file containing multiplexed tracings
	# [SAMPLESHEET.CSV] A 3 column csv (Lane,Sample,Index) formatted according to https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/bcl2fastq-direct
# To create SampleSheet.csv
#	
#	awk -F "," '{print "*",substr($2,1,length($2)-1),substr($3,1,length($3)-1)}' OFS="," samplesheet.csv | \sort | uniq 

BCL=${1}
SAMPLESHEET=${2}

if [[ ! -d ${BCL} ]]; then
echo "Please run this script in the directory that contains ${BCL}"
exit
fi

if [[ ! -f ${SAMPLESHEET} ]]; then
echo "Please run this script in the directory that contains ${SAMPLESHEET}"
exit
fi

# Mask bases needs to correspond to RunInfo.xml sheet in bcl file
#BASEMASK="Y151,I8N9,N8,Y141"
BASEMASK="Y151,I8,N8,Y151"


echo '#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --job-name=job_wrapper
#SBATCH --account=PAS1555


cd $SLURM_SUBMIT_DIR
pwd
source activate ~osu9900/anaconda2/envs/myenv2.7/


cellranger mkfastq --ignore-dual-index --id=mkfastq_'${SAMPLESHEET%.*}' --run=./'${BCL}' --use-bases-mask='${BASEMASK}' --csv='${SAMPLESHEET}'




' | sbatch



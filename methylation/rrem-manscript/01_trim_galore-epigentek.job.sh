#!/bin/bash


for FQR1 in *R1*.fastq.gz;
do

if [ -a "$FQR1" ];then


FQNAME=${FQR1%_R1*}

#CHECK FOR PREVIOUS TRIMMING:
if [ -a ${FQNAME}/${FQR1}_trimming_report.txt ];then
	printf ":::WARNING:::\n\tSKIPPING SAMPLE: ${FQNAME} !!!!\n"
else

printf "Submitting trim_galore for ${FQNAME}: "

echo '#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --job-name=trim_galore:'$FQNAME'
#SBATCH --output=job.trim_galore.'`date +%Y%m%d.%H%M%S`'.'$FQNAME'.log
#SBATCH --account=PASXXXX

date

cd $SLURM_SUBMIT_DIR


FQNAME='$FQNAME'
FQR1='$FQR1'
FQR2=$(ls ${FQNAME}*R2*)


if [ -d "./${FQNAME}/" ];then
	echo "output to directory trim_galore"
else
	mkdir ${FQNAME}
        echo "output to directory trim_galore"
fi


printf "Sample being processed is $FQNAME\n"
printf "\tRead 1 file: ${FQR1}\n"
printf "\tRead 2 file: ${FQR2}\n"

echo "Output directory is "$SLURM_SUBMIT_DIR"/${FQNAME}"

if [ -d ${FQNAME} ];then
	true
else
	mkdir ${FQNAME}
fi


/users/PAS0472/osu6369/tools/trim_galore_v0.4.0/trim_galore -q 20 --phred33 --fastqc \
	--gzip --rrbs --non_directional --paired --retain_unpaired -o ${SLURM_SUBMIT_DIR}/${FQNAME}/ ${FQR1} ${FQR2}



date
' | sbatch
fi

else
echo "No files exist matching "$FQR1
fi

done

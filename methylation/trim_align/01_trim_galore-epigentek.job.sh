#!/usr/bin/env bash


for FQR1 in *R1*.fastq.gz;
do

if [ -a "$FQR1" ];then


FQNAME=${FQR1%_R1*}

#CHECK FOR PREVIOUS TRIMMING:
if [ -a ${FQNAME}/${FQR1}_trimming_report.txt ];then
	printf ":::WARNING:::\n\tSKIPPING SAMPLE: ${FQNAME} !!!!\n"
else

printf "Submitting trim_galore for ${FQNAME}: "
group=$(echo $(groups) | cut -d" " -f1)


echo '
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4GB
#PBS -N trim_galore:'$FQNAME'
#PBS -S /bin/bash
#PBS -o job.trim_galore.'`date +%Y%m%d.%H%M%S`'.'$FQNAME'.log
#PBS -j oe
#PBS -W umask=0022
#PBS -A '${group}'

date

cd $PBS_O_WORKDIR


FQNAME='$FQNAME'
FQR1='$FQR1'
FQR2=$(ls ${FQNAME}*R2*)


#if [[ '$FQR1' == "*.gz" ]];then
#	zcat '$FQ' >$TMPDIR/'$FQNAME'
#	mkdir '$OUTDIR'
#else
#	cp '$FQ' $TMPDIR/'$FQNAME'
#	mkdir '$OUTDIR'
#fi

if [ -d "./${FQNAME}/" ];then
	echo "output to directory trim_galore"
else
	mkdir ${FQNAME}
        echo "output to directory trim_galore"
fi


printf "Sample being processed is $FQNAME\n"
printf "\tRead 1 file: ${FQR1}\n"
printf "\tRead 2 file: ${FQR2}\n"

echo "Output directory is "$PBS_O_WORKDIR"/${FQNAME}"

if [ -d ${FQNAME} ];then
	true
else
	mkdir ${FQNAME}
fi


/users/PAS0472/osu6369/tools/trim_galore_v0.4.0/trim_galore -q 20 --phred33 --fastqc --gzip --rrbs --non_directional --paired --retain_unpaired -o $PBS_O_WORKDIR/${FQNAME}/ ${FQR1} ${FQR2}



date
' | qsub
fi

else
echo "No files exist matching "$FQ
fi

done

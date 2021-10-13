echo '
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=10
#PBS -N job_wrapper
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=002
#PBS -A PAS1374


date

cd $PBS_O_WORKDIR
#source activate ~/anaconda2/envs/myenv2.7/
pwd
#python bam_dup.py -b '${1}'
module load samtools

echo '${1}'
#python subsample_30million/bam_dup_v2.py -b '${1}' -o pos

samtools collate '${1}' '${1%bam}'RANDOM
samtools view '${1%bam}'RANDOM.bam | head -n '${2}' > '${1%bam}'SUBSAMPLED.sam
samtools view -H '${1}' > '${1}'.tmp.header
cat '${1}'.tmp.header '${1%bam}'SUBSAMPLED.sam > '${1}'.tmp && mv '${1}'.tmp '${1%bam}'SUBSAMPLED.sam
\rm '${1}'.tmp.header '${1%bam}'RANDOM.bam
samtools view -b '${1%bam}'SUBSAMPLED.sam > ./subsample_pairwise_unalignedReads/'${1%bam}'SUBSAMPLED.'${2}'.bam
\rm '${1%bam}'SUBSAMPLED.sam
samtools sort ./subsample_pairwise_unalignedReads/'${1%bam}'SUBSAMPLED.'${2}'.bam > ./subsample_pairwise_unalignedReads/'${1%bam}'SUBSAMPLED.'${2}'.sorted.bam


#n=$(samtools view '${1}' | fgrep -f chrms.txt -w | wc -l)
#echo $n
#samtools index '${1}'
#samtools sort -o '${1%bam}'.sorted.bam '${1}'

#picard -Xms1g -Xmx4g CollectGcBiasMetrics I='${1}' O='${1%%_*}'_gc_metrics.txt CHART='${1%%_*}'_gc_metrics.pdf S='${1%%_*}'_gc_metrics.summary \
#	R=hg38.new.fa BS=true MAX_RECORDS_IN_RAM=100000 STOP_AFTER=$n
date
echo "Workflow complete"
' | qsub


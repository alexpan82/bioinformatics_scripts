echo '
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=16
#PBS -N job_wrapper
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=002
#PBS -A PAS1427


date

cd $PBS_O_WORKDIR

pwd
module load gnu/4.8.5
module load bedtools
ls '${1}' '${2}' 

ROINAME='${2#roi/}'
ROINAME=${ROINAME#refgene.}
ROINAME=${ROINAME%.hg38*}.bed

echo '${1}'.${ROINAME}.tmp
echo '${1}'.${ROINAME}.bed
echo '${1%.sorted.bam_CpG.txt}'.${ROINAME}

awk -F "\t" '"'"'{print $2 "\t" $3 "\t" $3}'"'"' '${1}' > '${1}'.${ROINAME}.tmp

paste '${1}'.${ROINAME}.tmp '${1}' | tail -n +2 > '${1}'.${ROINAME}.bed
\rm '${1}'.${ROINAME}.tmp

bedtools intersect -wo -bed -f 1E-9 -a '${1}'.${ROINAME}.bed -b '${2}' > '${1%.sorted.bam_CpG.txt}'.${ROINAME}

\rm '${1}'.${ROINAME}.bed

date
echo "Workflow complete"
' | qsub


#! /bin/bash
#takes an unsorted sam file (bam) with headers, sorts it into methylKit format, then runs r script to calculate the percent methylation. then runs python script to create file with only the percent methylation.

if [ $# -ne 1 ]
then

echo "Usage: genome_assembly"

else

ASSEMBLY=$1
for i in $(cat samples_into_methylkit.txt); do
#for i in test.bam; do
SAM=${i%%.*}


echo '
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=15
#PBS -l mem=60GB
#PBS -N methylkit_'${SAM%%.*bam}'
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=002
#PBS -A PAS1359

date
pwd

if [ -d $PBS_O_WORKDIR/mkitpermeth/ ];then
        echo "Outputting to "$PBS_O_WORKDIR/mkitpermeth/
else 
        echo "Making output directory "$PBS_O_WORKDIR/mkitpermeth/
        mkdir $PBS_O_WORKDIR/mkitpermeth/
fi

cd $PBS_O_WORKDIR

cp '${SAM}'.bam $TMPDIR

cd $TMPDIR

pwd
ls
echo '${ASSEMBLY}'



~osu5996/bin/samtools-0.1.10/samtools sort -m 1000000000 '${SAM}'.bam '${SAM}'.sort
date
pwd
ls

~osu5996/bin/samtools-0.1.10/samtools view -h -o '${SAM}'.sort.sam '${SAM}'.sort.bam

#cp '${SAM}'.sort.sam $PBS_O_WORKDIR/mkitpermeth/

#send sorted sam file to rscript to calculate percent methylation
echo "Creating percent methylation file for "'${SAM}'

module load R
echo "HI"
/usr/local/R/3.3.2/bin/R --slave --no-restore --file=/users/PAS0472/osu9900/RRBS_scripts/methylkit_analysis/methylKit2.R --args '${SAM}'.sort.sam '${SAM}' '${ASSEMBLY}'
#Rscript /fs/scratch/ccri0063/amal_cf/methylkit/methylKit2.R '${SAM}'.sort.sam '${SAM}' '${ASSEMBLY}'
echo "HI"

echo "Created percent methylation file as '${SAM}'_CpG.txt"
ls
cp '${SAM}'_CpG.txt $PBS_O_WORKDIR/mkitpermeth/'${SAM}'_CpG.txt

date
pwd

echo "Creating simpler percent methylation file for "'${SAM}'

python /fs/scratch/ccri0063/amal_cf/methylkit/permeth.methylKit.py '${SAM}'_CpG.txt '${SAM}'.permeth
ls
cp '${SAM}'.permeth $PBS_O_WORKDIR/mkitpermeth/

date
echo "Workflow complete"
' | qsub 
done

fi

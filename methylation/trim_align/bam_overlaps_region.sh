# 1st argument is a sorted bam file and 2nd argument is a bed file. 3rd is the outputname
# Intersects bam file locations with a region. Returns bed file that contain regions of the roi file that overlap with bam

TEMP=tmp_${3}.txt

echo '
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=12GB
#PBS -N bam_overlaps_region
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=002

date
cd $PBS_O_WORKDIR
cd ~/scratch/amal_cf/bismark_before_121817/sorted_bam
pwd
echo

module load gnu
module load bedtools
module load python


echo '${TEMP}'
bedtools intersect -wo -bed -abam '${1}' -b '${2}' > '${TEMP}'

awk '"'"'{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20}'"'"' '${TEMP}' > '${3}'

\rm '${TEMP}'

python unique_bamcounts_by_region.py '${3}'


echo "Workflow complete"

' | qsub

# 1st argument is a sorted bam file and 2nd argument is a bed file. 3rd is the outputname
#Reads_in_region.sh - counts unique reads that intersect at least one genomic range in an roi file
echo '
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=5
#PBS -N reads_in_region
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=002

date
cd $PBS_O_WORKDIR
pwd
echo

module load gnu
module load bedtools

bedtools intersect -abam '${1}' -b '${2}' -wa -u -f 1.0  > '${3}'

#Writes to text.txt
count=$(samtools view -c '${1}')
count_rpt=$(samtools view -c '${3}')

printf '${1%_R1*}' >> text1.txt
printf "\t" >> text1.txt
printf ${count} >> text1.txt
printf "\t" >> text1.txt
printf ${count_rpt} >> text1.txt
echo >> text1.txt
echo "Workflow complete"
' | qsub



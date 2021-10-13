#! /bin/bash
#takes an unsorted sam file (bam) with headers, sorts it into methylKit format, then runs r script to calculate the percent methylation. then runs python script to create file with only the percent methylation.

if [ $# -ne 1 ]
then

echo "Usage: genome_assembly"

else

ASSEMBLY=$1
for i in *bam; do
SAM=$(echo ${i} | sed 's/bam/sam/g')
BAM=${i}

echo '#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=methylkit_'${SAM%%.*bam}'
#SBATCH --account=PAS1829

date
pwd

if [ -d $SLURM_SUBMIT_DIR/mkitpermeth/ ];then
        echo "Outputting to "$SLURM_SUBMIT_DIR/mkitpermeth/
else 
        echo "Making output directory "$SLURM_SUBMIT_DIR/mkitpermeth/
        mkdir $SLURM_SUBMIT_DIR/mkitpermeth/
fi

cd $SLURM_SUBMIT_DIR

cp '${BAM}' $TMPDIR

cd $TMPDIR

pwd
ls
echo '${ASSEMBLY}'


module load samtools
samtools sort '${BAM}' > '${BAM%bam}'sorted.bam 
#date
#pwd
ls '${BAM%bam}'sorted.bam

#samtools view -h '${SAM}'.sorted.bam > '${SAM}'.sort.sam 
samtools view -h '${BAM%bam}'sorted.bam  > '${SAM}'.sort.sam

#send sorted sam file to rscript to calculate percent methylation
echo "Creating percent methylation file for "'${BAM}'

source activate ~osu9900/anaconda2/envs/r-methylkit/

Rscript ~osu9900/RRBS_scripts/methylkit_analysis/methylKit2.R '${SAM}'.sort.sam '${BAM}' '${ASSEMBLY}'

echo "Created percent methylation file as '${BAM}'_CpG.txt"
ls
cp '${BAM}'_CpG.txt $SLURM_SUBMIT_DIR/mkitpermeth/'${BAM}'_CpG.txt

date
pwd


date
echo "Workflow complete"
' | sbatch 
done

fi

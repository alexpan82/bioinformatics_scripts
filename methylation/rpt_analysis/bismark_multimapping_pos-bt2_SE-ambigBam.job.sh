#!/bin/bash
# This script takes 3 arguements. 1st is the genome name (see below) and 2nd is a file containing the files you wish to align (each on a new line), and 3rd is the upper-limit on the number of reads you'd like to return for each multimapped read

#The human genome that will be referenced is taken in as an argument
if [[ "$1" == "hg19" ]]; then
    GENOME=$1
    BIS_PATH=/nfs/06/osu6369/PAS0472/proj11/data/reference/$GENOME/bismark_bt2_index/

elif [[ "$1" == "hg18" ]]; then
    GENOME=$1
    BIS_PATH=/nfs/06/osu6369/PAS0472/proj11/data/reference/$GENOME/bismark_bt2_index/

elif [[ "$1" == "mm9" ]]; then
    GENOME=$1
    BIS_PATH=/nfs/06/osu6369/PAS0472/proj11/data/reference/$GENOME/bismark_bt2_index/

elif [[ "$1" == "mm10" ]]; then
    GENOME=$1
    BIS_PATH=/users/PAS0472/osu6369/PAS0472/osu7905/references/genomes/mouse/mm10/bismark_bt2_index/
elif [[ "$1" == "GRCh38" ]]; then
    GENOME=$1
    #latest bismark version v0.17 built GRCh38 genome WITHOUT unplaced chr
    BIS_PATH=/users/PAS0472/osu6369/PAS0472/osu7905/references/genomes/human/GRCh38/bismark_bt2_index-current/
    #latest bismark version v0.17 built full GRCh38 genome
    #BIS_PATH=/fs/project/PAS0472/osu7905/references/genomes/human/GRCh38/bismark_bt2_index-current/
    #temp lustre location
    #BIS_PATH=/fs/scratch/osu7905/references/genomes/human/GRCh38/bismark_bt2_index/
elif [[ "$1" == "bc_AU1054" ]]; then
    GENOME=$1
    BIS_PATH=/fs/scratch/ccri0063/amal_cf/reference/bc_AU1054/bsgenome2/
else
       echo "Input library used to align (hg18, hg19, GRCh38, mm9, or mm10)";
fi



# test directory structure
if [[ ! -z "$1" ]]; then


    if [ ! -d ../bismark ]; then
	mkdir ../bismark
    fi

    DIR=`pwd`;


#loop through read1(R1).fq.gz files
for FASTQ in $(cat ${2}); do


group=$(echo $(groups) | cut -d" " -f1)

printf "Alignment job for ${FASTQ} is: "
echo '
#PBS -N job.'$FASTQ'.'$GENOME'.bismark-bt2_PE.job
#PBS -l nodes=1:ppn=28
#PBS -l walltime=10:00:00
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=022
#PBS -A PAS1359


date
echo "Reference Genome used for alignment is '$GENOME'"
echo "Creating Bismark Alignment"
cd $PBS_O_WORKDIR

#import genome path
BIS_PATH='${BIS_PATH}'

#import fastq basename

printf "Processing SE files with basename: '${FASTQ}'\n"

#GET FASTQ

FASTQGZ='${FASTQ}'
FASTQ=${FASTQGZ%.gz}
zcat ${FASTQGZ} > $TMPDIR/${FASTQ}
printf "Single-end file: $(ls ${TMPDIR}/${FASTQ})\n"


#MOVE TO TMPDIR, COUNT READS, ALIGN
cd $TMPDIR
printf "Total non-duplicated reads that pass filter: \n"
n1=($( wc -l ${FASTQ}));
n1=$(( n1/4 ))
echo ${n1}


###
#   :::TMP FIX::: Use a dump directory to dump all output of bismark!
###

#make dump dir:
printf "bismark dump dir is: ${PBS_O_WORKDIR}/../bismark\n"
mkdir -p ${PBS_O_WORKDIR}/../bismark

module load samtools

# Options for bismark are changed here, currently set up for paired end alignment, see manual for other options
#latest version of bismark

echo -e "Upper-limit of the # of locations returned for each multimapped read:\t"
printf '${3}'

echo ":::Bismark Parameters::: /users/PAS0472/osu9861/tools/bismark_multimapping_locations/bismarkMultimappingPositions --multicore 4 -N 1 -p 2 --non_directional --score_min L,-0.6,-0.6 -o ${PBS_O_WORKDIR}/../bismark --genome ${BIS_PATH} ${FASTQ} --limit 20"

#/users/PAS0472/osu9861/tools/bismark_multimapping_locations/bismarkMultimappingPositions --bowtie2 --path_to_bowtie ~osu6369/tools/bowtie2-current/ --multicore 4 --unmapped --ambiguous --ambig_bam --non_directional -N 1 -p 2 --score_min L,-0.6,-0.6 -o ${PBS_O_WORKDIR}/../bismark --limit 20 ${BIS_PATH} ${FASTQ}

/users/PAS0472/osu9861/tools/bismark_multimapping_locations/bismarkMultimappingPositions --bowtie2 --path_to_bowtie ~osu6369/tools/bowtie2-current/ --multicore 4 -N 1 -p 2 --non_directional --score_min L,-0.6,-0.6 -o ${PBS_O_WORKDIR}/../bismark --genome ${BIS_PATH} ${FASTQ} --limit_multimapping '${3}'


printf "Finished with Bismark pipeline on '$FASTQ'\n"


' |qsub

   done;

else
    echo "Nothing submitted.";
fi

date

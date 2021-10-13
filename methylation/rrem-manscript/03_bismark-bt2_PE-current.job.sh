#!/bin/bash


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
for FASTQ in *R1*val*.fq.gz;do
	if [ -a $FASTQ ];then
		FQNAME=${FASTQ%_R1*}
	else
		echo 'No files matching '$FASTQ' found'
	fi


group=$(echo $(groups) | cut -d" " -f1)


printf "Alignment job for ${FQNAME} is: "
echo '#!/bin/bash
#SBATCH --time=25:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=25
#SBATCH --job-name=job.'$FQNAME'.'$GENOME'.bismark-bt2_PE.job
#SBATCH --output=job.'$FQNAME'.'$GENOME'.bismark-bt2_PE.'`date +%Y%m%d.%H%M%S`'.'$FQNAME'.log
#SBATCH --account=PASXXXX


date
echo "Reference Genome used for alignment is '$GENOME'"
echo "Creating Bismark Alignment"
cd $SLURM_SUBMIT_DIR

#import genome path
BIS_PATH='${BIS_PATH}'

#import fastq basename
FQNAME='${FQNAME}'
printf "Processing PE files with basename: ${FQNAME}\n"

#GET R1 and R2 files
FASTQR1GZ=$(ls '$FQNAME'*R1*.fq.gz )
FASTQR1=${FASTQR1GZ%.gz}
zcat ${FASTQR1GZ} > $TMPDIR/${FASTQR1}
printf "Read 1 file: $(ls ${TMPDIR}/${FASTQR1})\n"

FASTQR2GZ=$( ls '$FQNAME'*R2*.fq.gz )
FASTQR2=${FASTQR2GZ%.gz}
zcat ${FASTQR2GZ} > $TMPDIR/${FASTQR2}
printf "Read 2 file: $(ls ${TMPDIR}/${FASTQR2})\n"

#MOVE TO TMPDIR, COUNT READS, ALIGN
cd $TMPDIR
printf "Total non-duplicated reads that pass filter: \n"
n1=($( wc -l ${FASTQR1}));
n1=$(( n1/4 ))
n2=($( wc -l ${FASTQR2}));
n2=$(( n2/4 ))
echo ${n1}
echo ${n2}

###
#   :::TMP FIX::: Use a dump directory to dump all output of bismark!
###

#make dump dir:
printf "bismark dump dir is: ${SLURM_SUBMIT_DIR}/../bismark\n"
mkdir -p ${SLURM_SUBMIT_DIR}/../bismark

module load samtools

# Options for bismark are changed here, currently set up for paired end alignment, see manual for other options
#~osu6369/tools/bismark-current/bismark --bowtie2 --path_to_bowtie ~osu6369/tools/bowtie2-current/ --multicore 1 --unmapped --ambiguous --non_directional -X 1000 -N 1 -p 2 --score_min L,-0.6,-0.6 -o ${SLURM_SUBMIT_DIR}/../bismark ${BIS_PATH} -1 ${FASTQR1} -2 ${FASTQR2}
#latest version of bismark
echo ":::Bismark Parameters::: /users/PAS0472/osu6369/tools/bismark_v0.17.0/bismark --bowtie2 --path_to_bowtie ~osu6369/tools/bowtie2-current/ --multicore 1 --non_directional -X 1000 -N 1 p 2 --score_min L,-0.6,-0.6 -o ${SLURM_SUBMIT_DIR}/../bismark ${BIS_PATH} -1 ${FASTQR1} -2 ${FASTQR2}"
#/users/PAS0472/osu6369/tools/bismark_v0.17.0/bismark --bowtie2 --path_to_bowtie ~osu6369/tools/bowtie2-current/ --multicore 1 --unmapped --ambiguous --non_directional -X 1000 -N 1 -p 1 --score_min L,-0.6,-0.6 -o ${SLURM_SUBMIT_DIR}/../bismark ${BIS_PATH} -1 ${FASTQR1} -2 ${FASTQR2}
/users/PAS0472/osu6369/tools/bismark_v0.17.0/bismark --bowtie2 \
    --path_to_bowtie ~osu6369/tools/bowtie2-current/ --multicore 4 \
    --non_directional -X 1000 -N 1 -p 4 --score_min L,-0.6,-0.6 \
    -o ${SLURM_SUBMIT_DIR}/../bismark ${BIS_PATH} -1 ${FASTQR1} -2 ${FASTQR2}


#printf "Creating methylation calls using \"bismark_methylation_extractor\" function:\n"
#make directory
#if [ -d ${SLURM_SUBMIT_DIR}/../methExt ];then
#    true
#else
#    mkdir ${SLURM_SUBMIT_DIR}/../methExt
#fi

# Full methylation extraction:
#printf "Performing "bismark_methylation_extractor" using PE aligned file: \n\t $(ls ${SLURM_SUBMIT_DIR}/../bismark/${FQNAME}*bismark_bt2_pe.sam)\n"
#old bismark version
#~osu6369/tools/bismark-current/bismark_methylation_extractor -p --bedGraph --counts --genome_folder ${BIS_PATH} --comprehensive -o ${SLURM_SUBMIT_DIR}/../methExt ${SLURM_SUBMIT_DIR}/../bismark/${FQNAME}*bismark_bt2_pe.bam
#/users/PAS0472/osu6369/tools/bismark_v0.17.0/bismark_methylation_extractor -p --bedGraph --no_overlap --counts --gzip --genome_folder ${BIS_PATH} --comprehensive -o ${SLURM_SUBMIT_DIR}/../methExt ${SLURM_SUBMIT_DIR}/../bismark/${FQNAME}*bismark_bt2_pe.bam


printf "Finished with Bismark pipeline on '$FQNAME'\n"

#OLD FILE MIGRATION
#echo "Moving files to directories..."

#gzip '$FQNAME'_bismark.sam

#cd $SLURM_SUBMIT_DIR;
#cp $TMPDIR/'$FASTQ'_bismark.sam.gz ../bismark/sam
#cp $TMPDIR/*.txt ../bismark/meth_report
#cp $TMPDIR/*.bam ../bismark/bam
#echo "Done moving files."

' | sbatch

   done;

else
    echo "Nothing submitted.";
fi

date

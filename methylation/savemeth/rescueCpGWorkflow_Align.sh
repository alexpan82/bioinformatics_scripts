### This script is a workflow to resuce methylation from multimapped and unaligned reads
### Run this script in a directory containing all fastq files
### Run rescueCpGWorkflow_Save.sh after this script

### Workflow ###
## 1 Combine multimapped & unaligned PE fastqs into SE fastqs (Please use trimmed fastqs)
	# Outputted from PE Bismark alignment (*ambiguous.fq.gz & *unmapped.fq.gz)
## 2 Obtain ambiguous read locations
	# Run BismarkMultimappingPositions on combined fastqs
## 3 Run Bismark SE alignment on the files curcurrently
## 4 Save CpGs and return a tsv file (#chr \t CpGposition \t coverage \t number of Cs \t methylation ratio (numC/cov))

### Script takes sample name prefixes in a text file, with each name being on a new line
	### For example, the fastq file test1_V1B_comb_SE_reads.ambiguous.fq.gz should be denoted as 'test1' in the text file
	### for i in *fq.gz; do echo ${i%%_*}; done | \sort | uniq > samples.txt

### Usage: bash rescueCpGWorkflow_Align.sh [sample list] [path to bismark indexed & converted genome to align to] [limit on multimap locations (int)]
### Ex: bash rescueCpGWorkflow_Align.sh samples.txt GRCh38 20
	# This finds out to 40 multimapping locations for each mate

if [[ ${1} == "" ]]; then
	echo "Please supply a text file containing sample name identifiers on each new line"
	echo "Usage: bash rescueCpGWorkflow_Align.sh [sample list] [path to bismark indexed & converted genome to align to] [limit on multimap locations (int)]"
	exit
fi

#The human genome that will be referenced is taken in as an argument
if [[ "$2" == "hg19" ]]; then
	GENOME=$2
	BIS_PATH=/nfs/06/osu6369/PAS0472/proj11/data/reference/$GENOME/bismark_bt2_index/

elif [[ "$2" == "hg18" ]]; then
	GENOME=$2
	BIS_PATH=/nfs/06/osu6369/PAS0472/proj11/data/reference/$GENOME/bismark_bt2_index/

elif [[ "$2" == "mm9" ]]; then
	GENOME=$2
	BIS_PATH=/nfs/06/osu6369/PAS0472/proj11/data/reference/$GENOME/bismark_bt2_index/

elif [[ "$2" == "mm10" ]]; then
	GENOME=$2
	BIS_PATH=/users/PAS0472/osu6369/PAS0472/osu7905/references/genomes/mouse/mm10/bismark_bt2_index/
elif [[ "$2" == "GRCh38" ]]; then
	GENOME=$2
	#latest bismark version v0.17 built GRCh38 genome WITHOUT unplaced chr
	BIS_PATH=/users/PAS0472/osu6369/PAS0472/osu7905/references/genomes/human/GRCh38/bismark_bt2_index-current/
elif [[ "$2" == "bc_AU1054" ]]; then
	GENOME=$2
	BIS_PATH=/fs/scratch/ccri0063/amal_cf/reference/bc_AU1054/bsgenome2/
else
	echo "Input library used to align (hg18, hg19, GRCh38, mm9, or mm10)";
fi

if [[ ${3} == "" ]]; then
	echo "Please supply an integer to specify the top number of multimapping locations for each read to report (/2)"
	echo "Usage: bash findConsensusCpG.sh [sample list] [genome to align to] [limit on multimap locations (int)]"
	exit
fi

# test directory structure
if [ ! -d ../bismark_SEalign ]; then
	mkdir -p ../bismark_SEalign
fi

if [ ! -d ../bismark_multimapLOC ]; then
	mkdir -p ../bismark_multimapLOC
fi


# Run job for every sample
for NAME in $(cat $1); do

	# Names of trimmed mates
	# a for ambiguous, u for unmapped/unaligned
	aR1NAME=$(ls ${NAME}*R1*ambiguous*fq.gz)
	aR2NAME=$(ls ${NAME}*R2*ambiguous*fq.gz)
	uR1NAME=$(ls ${NAME}*R1*unmapped*fq.gz)
	uR2NAME=$(ls ${NAME}*R2*unmapped*fq.gz)
	

	echo '
	#PBS -l walltime=20:00:00
	#PBS -l nodes=1:ppn=28
	#PBS -l mem=112GB
	#PBS -N savemeth_bmp_'${NAME}'
	#PBS -S /bin/bash
	#PBS -j oe
	#PBS -W umask=002
	#PBS -A PAS1359

	date

	# This is fastq/
	cd $PBS_O_WORKDIR
	pwd

	# This may need to be commented out
	#source activate ~/anaconda2/envs/myenvironment/
	#module load gnu/4.8.5
	#module load bedtools
	#module load samtools

	aR1NAME=$(ls '${aR1NAME}')
	aR2NAME=$(ls '${aR2NAME}')
	uR1NAME=$(ls '${uR1NAME}')
	uR2NAME=$(ls '${uR2NAME}')


	################################################### Combine fastqs

	COMBFASTQ=${aR1NAME%%_*}_V1B_comb_SE_reads.NOTuniquemap.fq.gz
	
	echo "Combining" ${aR1NAME} ${aR2NAME} ${uR1NAME} ${uR2NAME}
	cat ${aR1NAME} ${aR2NAME} ${uR1NAME} ${uR2NAME} > $COMBFASTQ
	
	
	################################################### Obtain ambiguous read locations and perform Bismark SE alignment in a separate job

	FASTQGZ=$COMBFASTQ
	FASTQ=${FASTQGZ%.gz}

	echo "Reference Genome used for alignment is '$GENOME'"
		#import genome path
	BIS_PATH='${BIS_PATH}'
	printf "Processing SE files with basename: ${COMBFASTQ} \n"

	zcat ${FASTQGZ} > $TMPDIR/${FASTQ}

	printf "Total non-duplicated reads that pass filter: \n"
	n1=($( wc -l $TMPDIR/${FASTQ}));
	n1=$(( n1/4 ))
	echo ${n1}
	\rm $TMPDIR/${FASTQ}

	# Options for bismark are changed here, currently set up for paired end alignment, see manual for other options
	#latest version of bismark

	echo -e "Upper-limit of the # of locations returned for each multimapped read:\t" '${3}'

	echo

	echo "Bismark SE Alignment Parameters::: /users/PAS0472/osu6369/tools/bismark_v0.17.0/bismark \
	--bowtie2 --path_to_bowtie ~osu6369/tools/bowtie2-current/ --unmapped --ambiguous --non_directional \
	-N 1 -p 2 --score_min L,-0.6,-0.6 -o ${PBS_O_WORKDIR}/../bismark_SEalign ${BIS_PATH} ${FASTQGZ}"

	echo "
		#PBS -l walltime=25:00:00
		#PBS -l nodes=1:ppn=28
		#PBS -l mem=112GB
		#PBS -N savemeth_SEalign_'${NAME}'
		#PBS -S /bin/bash
		#PBS -j oe
		#PBS -W umask=002
		#PBS -A PAS1359

		cd $PBS_O_WORKDIR
		pwd
		echo ${BIS_PATH} ${FASTQGZ}
		/users/PAS0472/osu6369/tools/bismark_v0.17.0/bismark --bowtie2 --path_to_bowtie ~osu6369/tools/bowtie2-current/ \
		--unmapped --ambiguous --non_directional -N 1 -p 2 --score_min L,-0.6,-0.6 -o ${PBS_O_WORKDIR}/../bismark_SEalign ${BIS_PATH} ${FASTQGZ}

		################################################### Retain the fastq names of multimapped and unaligned reads

		zcat ${aR1NAME} ${aR2NAME} | fgrep "@" > ${PBS_O_WORKDIR}/../bismark_SEalign/'${NAME}'_ambig_names.txt
		zcat ${uR1NAME} ${uR2NAME} | fgrep "@" > ${PBS_O_WORKDIR}/../bismark_SEalign/'${NAME}'_unaligned_names.txt
		sed -i '"'"'s/@//g'"'"' ${PBS_O_WORKDIR}/../bismark_SEalign/'${NAME}'_ambig_names.txt
		sed -i '"'"'s/@//g'"'"' ${PBS_O_WORKDIR}/../bismark_SEalign/'${NAME}'_unaligned_names.txt

		" | qsub

	echo 

	echo ":::BismarkMultimappingLocations Parameters::: /users/PAS0472/osu9861/tools/bismark_multimapping_locations/bismarkMultimappingPositions \
		--un --ambiguous --bowtie2--path_to_bowtie ~osu6369/tools/bowtie2-current/ --multicore 4 -N 1 -p 2 --non_directional \
			--score_min L,-0.6,-0.6 -o ${PBS_O_WORKDIR}/../bismark_multimapLOC --genome ${BIS_PATH} ${FASTQGZ} --limit_multimapping" '${3}'

	/users/PAS0472/osu9861/tools/bismark_multimapping_locations/bismarkMultimappingPositions --un --ambiguous --bowtie2 \
		--path_to_bowtie ~osu6369/tools/bowtie2-current/ --multicore 4 -N 1 -p 2 --non_directional \
			--score_min L,-0.6,-0.6 -o ${PBS_O_WORKDIR}/../bismark_multimapLOC --genome ${BIS_PATH} ${FASTQGZ} --limit_multimapping '${3}'

	cd ${PBS_O_WORKDIR}/../bismark_multimapLOC
	\rm ${COMBFASTQ%%_*}*NOTuniquemap_bismark_bt2.bam ${COMBFASTQ%%_*}*NOTuniquemap_bismark_bt2_SE_report.txt
	\rm ${COMBFASTQ%%_*}*NOTuniquemap.fq.gz_ambiguous_reads.fq.gz ${COMBFASTQ%%_*}*NOTuniquemap.fq.gz_unmapped_reads.fq.gz 

  	printf "Finished with Bismark pipeline on ${COMBFASTQ} \n"


 

	echo "Workflow complete"

	' | qsub

done

# This script combines all fastq alliqouts with the SAME LIBRARY # into 1 file or merges bam files together using samtools
# Assumes paired-end gzipped fastqs
# Ex: bash comb_alliquots.sh fastq 
# Ex: bash comb_alliquots.sh bam
# Ex: bash comb_alliquots.sh weldr
	# This combines trim* align* in the output folder into 1 directory

if [[ ${1} == "" ]]
then 
	echo "Please input either fastq or bam as argument."
	exit


elif [[ ${1} == "fastq" ]]
then
	uniqnames=$(for name in *fastq.gz; do
		#echo ${name%%_*}
		echo ${name%_[1-5]_S*}
	done | \sort | uniq)

	for name in $uniqnames; do
		if [ $(ls ${name}* | wc -l) -gt 2 ]
		then
			echo 'Combining' $(ls ${name}*R1* | \sort) '>' ${name}_comb$(($(ls ${name}* | wc -l)/2))_R1.fastq.gz
			cat $(ls ${name}*R1* | \sort) > ${name}_comb$(($(ls ${name}* | wc -l)/2))_R1.fastq.gz
			echo 'Combining' $(ls ${name}*R2* | \sort) '>' ${name}_comb$(($(ls ${name}* | wc -l)/2))_R2.fastq.gz
			cat $(ls ${name}*R2* | \sort) > ${name}_comb$(($(ls ${name}* | wc -l)/2))_R2.fastq.gz
		
		fi
	done

	echo 'Alliquot combining completed successfully'

elif [[ ${1} == "bam" ]]
then
	echo 'Need to been updated. Functionality does not exist yet'

elif [[ ${1} == "weldr" ]]
then
	module load samtools

	uniqnames=$(for name in $(ls -d */); do
		echo ${name%_[1-5]_S*}
	done | \sort | uniq)	

	for name in $uniqnames; do
		if [ $(ls -d ${name}* | wc -l) -gt 1 ]
		then
			newdirname=${name}_comb$(ls -d ${name}* | wc -l)
			mkdir $newdirname
			echo 'Combining' $(ls ${name}*/trim_R1.fastq.gz | \sort) '>' ${newdirname}/trim_R1.fastq.gz
			cat $(ls ${name}*/trim_R1.fastq.gz | \sort) > ${newdirname}/trim_R1.fastq.gz
			echo 'Combining' $(ls ${name}*/trim_R2.fastq.gz | \sort) '>'${newdirname}/trim_R2.fastq.gz
			cat $(ls ${name}*/trim_R2.fastq.gz | \sort) > ${newdirname}/trim_R2.fastq.gz
			echo 'Combining' $(ls ${name}*/align.bam | \sort) '>' ${newdirname}/align.bam
			samtools merge ${newdirname}/align.bam $(ls ${name}*/align.bam | \sort)
			
		fi
	done

else
	echo "Please input either fastq, bam, or weldr as argument."
	exit
fi



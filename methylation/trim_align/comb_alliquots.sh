# This script combines all fastq alliqouts into 1 file or merges bam files together using samtools
# Assumes paired-end gzipped fastqs
# Ex: bash comb_alliquots.sh fastq 
# Ex: bash comb_alliquots.sh bam

if [[ ${1} == "" ]]
then 
	echo "Please input either fastq or bam as argument."
	exit


elif [[ ${1} == "fastq" ]]
then
	uniqnames=$(for name in *fastq.gz; do
		echo ${name%_[1-5]_S*}
	done | \sort | uniq)

	for name in $uniqnames; do
		if [ $(ls ${name}* | wc -l) -gt 2 ]
		then
			echo 'Combining' $(ls ${name}*R1* | \sort) '>' ${name}_comb$(($(ls ${name}* | wc -l)/2))_SC_L00C_R1.fastq.gz
			cat $(ls ${name}*R1* | \sort) > ${name}_comb$(($(ls ${name}* | wc -l)/2))_SC_L00C_R1.fastq.gz
			echo 'Combining' $(ls ${name}*R2* | \sort) '>' ${name}_comb$(($(ls ${name}* | wc -l)/2))_SC_L00C_R2.fastq.gz
			cat $(ls ${name}*R2* | \sort) > ${name}_comb$(($(ls ${name}* | wc -l)/2))_SC_L00C_R2.fastq.gz
		
		fi
	done

	echo 'Alliquot combining completed successfully'

elif [[ ${1} == "bam" ]]
then
	echo 'Need to been updated. Functionality does not exist yet'

else
	echo "Please input either fastq or bam as argument."
	exit
fi



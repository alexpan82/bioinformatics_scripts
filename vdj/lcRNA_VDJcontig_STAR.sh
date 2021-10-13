# Builds VDJ contigs for reads with full CDR3
# Takes read names, trimmed PE reads, and output name
# Ex: bash lcRNA_VDJcontig.sh db-pass.tsv trim_R1.fastq.gz trim_R2.fastq.gz output

if [ ${1} != "-" ]; then
	NAMES=$(readlink -f ${1})
else
	NAMES=${1}
fi
R1=$(readlink -f ${2})
R2=$(readlink -f ${3})
OUTNAME=${4}

# Make output directory
#mkdir ${OUTNAME}
#cd ${OUTNAME}

# Get reads from trimmed files
#python ~osu9900/techManScripts/fastq_names.py ${R1} ${R2} ${NAMES}

# Get reads from trimmed files
if [ NAMES != "-" ]; then
	python ~osu9900/techManScripts/fastq_names.py ${R1} ${R2} ${NAMES}
	R1RUN="subset_R1.fastq"
	R2RUN="subset_R2.fastq"	
else
	R1RUN=${R1}
	R2RUN=${R2}	
fi

# Assemble
#/users/PAS0472/osu8725/tools/trinityrnaseq-Trinity-v2.8.2/Trinity --seqType fq --max_memory 50G \
#	--left ${R1RUN} --right ${R2RUN} --CPU 12

#seqtk seq -r trinity_out_dir/Trinity.fasta > trinity_out_dir/Trinity_revcomp.fasta

# Re-align reads to contigs
	# Prepare index

# Run igblast
igblastdir=~osu8725/tools/igblast/
IGOUT="igblast_output/Trinity"
if [ ! -d ${IGOUT%%/*} ]; then mkdir ${IGOUT%%/*}; fi

~osu8725/tools/igblast/bin/igblastn -query trinity_out_dir/Trinity.fasta \
	-num_threads 40 \
	-domain_system imgt \
	-ig_seqtype Ig \
	-organism human \
	-auxiliary_data ${igblastdir}/optional_file/human_gl.aux \
	-outfmt '7 std qseq sseq btop' \
	-out ${IGOUT}.fmt7

~osu8725/anaconda3/bin/MakeDb.py igblast -i ${IGOUT}.fmt7 \
	-s trinity_out_dir/Trinity.fasta \
	-r ~osu8725/tools/igblast/database/IMGT/IG[VDJ].fasta \
	--regions --scores --failed
for i in ${IGOUT}*.tab; do mv $i ${i%.tab}.tsv; done

# Figure out which assembly file has the CDR3 in the top strand orientation
CDR3=$(cut -f 29 ${IGOUT}*pass.tsv | tail -n +2 | \sort | uniq -c | \sort -k1nr | head -n 1 | awk '{print $2}')

if [ $(grep ${CDR3} trinity_out_dir/Trinity.fasta | wc -l) -eq 0 ]; then
	REF="./trinity_out_dir/Trinity_revcomp.fasta"
else
	REF="./trinity_out_dir/Trinity.fasta"
fi

echo $REF

# Create index
mkdir index-STAR

~osu8725/tools/STAR-2.6.0a/source/STAR --runMode genomeGenerate --genomeSAindexNbases 1 --runThreadN 5 --genomeDir ./index-STAR \
	--genomeFastaFiles ${REF}

# Align
mkdir align-STAR
STAROUT="align-STAR"

~osu8725/tools/STAR-2.6.0a/source/STAR --runThreadN 5 --genomeSAindexNbases 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignIntronMin 2 \
	--genomeDir ./index-STAR/ \
	--outSAMunmapped Within --readFilesCommand cat \
	--readFilesIn ${R1RUN} ${R2RUN} \
	--outFileNamePrefix ${STAROUT}/ --outSAMtype BAM SortedByCoordinate
~osu8725/tools/toolbelt/src/samtools/samtools index ./${STAROUT}/Aligned.sortedByCoord.out.bam
~osu8725/tools/toolbelt/src/samtools/samtools mpileup -f ./${REF} ./${STAROUT}/Aligned.sortedByCoord.out.bam > ./${STAROUT}/mpileup.txt
~osu8725/tools/toolbelt/src/samtools/samtools mpileup -Q 0 -f ./${REF} ./${STAROUT}/Aligned.sortedByCoord.out.bam > ./${STAROUT}/mpileup_raw.txt
cat ./${STAROUT}/mpileup_raw.txt | python /fs/project/PAS0472/toolbelt/scripts/parse_mpileup.py > ./${STAROUT}/mpileup_raw_cov.txt




cd ..
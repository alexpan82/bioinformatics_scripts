# Builds VDJ contigs for reads with full CDR3 provided in mixcr clones.txt file
# Takes read names, trimmed PE reads, and output name
# Ex: bash neb_mixcr_trinity.sh mixcr_align.txt trim_R1.fastq.gz trim_R2.fastq.gz mixcr_clones.txt

#ALIGN=$(readlink -f ${1})
#R1=$(readlink -f ${2})
#R2=$(readlink -f ${3})
ALIGN=$(ls ${1})
R1=$(ls ${2})
R2=$(ls ${3})
CLONES=${4}
CHAIN='IGHV'
NUM="10"

# Get top NUM clones from CHAIN
echo "Get top NUM clones from CHAIN"
CDR3LIST=$(echo -e "top${NUM}_${CHAIN}_CDR3s.txt")
grep ${CHAIN} ${CLONES} | head -n ${NUM} | cut -f 4 > ${CDR3LIST}
grep -w -f ${CDR3LIST} ${ALIGN} | awk -F "\t" '{print $NF}' | awk '{print $1}' >> names.txt

# Get reads from trimmed files
echo "Get reads from trimmed files"
python /users/PAS0472/osu9900/techManScripts/fastq_names.py ${R1} ${R2} names.txt
R1RUN="subset_R1.fastq"
R2RUN="subset_R2.fastq"	

# Assemble
/users/PAS0472/osu8725/tools/trinityrnaseq-Trinity-v2.8.2/Trinity --seqType fq --max_memory 50G \
	--left ${R1RUN} --right ${R2RUN} --CPU 12

seqtk seq -r trinity_out_dir/Trinity.fasta > trinity_out_dir/Trinity_revcomp.fasta

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

if [ $(grep ${CDR3} trinity_out_dir/Trinity.fasta | wc -l) -eq 0 ]; then
	REF="./trinity_out_dir/Trinity_revcomp.fasta"
else
	REF="./trinity_out_dir/Trinity.fasta"
fi

echo $REF
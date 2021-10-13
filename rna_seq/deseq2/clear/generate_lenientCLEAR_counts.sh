# Generates all samples needed for lenient CLEAR comparisons
# Assumes that counts tables and CLEAR passed genes are in above directory
# Ie: ../counts.rep.tsv and ../*/passed_genes.txt and ../coldata.txt

# Argurments:
# 1) The column header of the ../coldata.txt that separates the 2 groups
# 2) Name of group 1
# 3) Name of group 2
# 4) Optional: Additional word to grep for while subsetting samples
#              Can also be another samples.txt file that explicitly tells which samples to subset
#              txt file must have each sample on a new line. Has the same structure as samples.txt

CATEGORY=${1}
GROUP1=${2}
GROUP2=${3}
CATEGORY2=${4}

if [[ ${CATEGORY2} = "" ]]; then
	OUTNAME=${GROUP1}_vs_${GROUP2}
else
	OUTNAME=${GROUP1}_vs_${GROUP2}_${CATEGORY2%.*}
fi

if [ ! -d ${OUTNAME} ]; then
	mkdir ${OUTNAME}
fi

##### Step 1: Create subsampled coldata.txt and samples.txt file
echo "Creating subsampled coldata.txt and samples.txt file"

CAT_NCOL=$(for i in $(head -n 1 ../coldata.txt); do echo $i; done | \
	grep -n -w ${CATEGORY} | awk -F ":" '{print $1}')

head -n 1 ../coldata.txt > ${OUTNAME}/coldata.txt

if [ -f "${CATEGORY2}" ]; then
	awk -v ncol=${CAT_NCOL} -v grpa=${GROUP1} -v grpb=${GROUP2} \
		'{if($(ncol) == grpa || $(ncol) == grpb){print $0}}' ../coldata.txt | \
		grep -w -f ${CATEGORY2} >> ${OUTNAME}/coldata.txt

elif [ -z ${CATEGORY2} ]; then
	awk -v ncol=${CAT_NCOL} -v grpa=${GROUP1} -v grpb=${GROUP2} \
		'{if($(ncol) == grpa || $(ncol) == grpb){print $0}}' ../coldata.txt >> ${OUTNAME}/coldata.txt

else
	awk -v ncol=${CAT_NCOL} -v grpa=${GROUP1} -v grpb=${GROUP2} \
		'{if($(ncol) == grpa || $(ncol) == grpb){print $0}}' ../coldata.txt | \
		grep -w ${CATEGORY2} >> ${OUTNAME}/coldata.txt
fi
cut -f 1 ${OUTNAME}/coldata.txt | tail -n +2 > ${OUTNAME}/samples.txt

##### Step 2: Subsample counts.rep.tsv
echo "Subsampling counts.rep.tsv"

bash /users/PAS0472/osu9900/clonality_scripts/get_tsv_columns.sh \
	../counts.tsv ${OUTNAME}/samples.txt ${OUTNAME}/counts.tsv

##### Step 3: Run lenient CLEAR
echo "Running lenient CLEAR"
python /users/PAS0472/osu9900/rna_scripts/deseq2/clear/CLEAR_bounded_expression.py \
	lenient 0 ${GROUP1} ${GROUP2} ${CATEGORY} ${OUTNAME}/counts.rep.tsv ${OUTNAME}/coldata.txt \
	> ${OUTNAME}/clear.counts.tsv

echo "DONE. Please see ${OUTNAME} for all outputted files"

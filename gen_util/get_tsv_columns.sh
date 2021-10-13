# Given a text file with names of column headers in counts table, will return those headers in the order given
	# Will automatically reutrn the 1st column (gene names) so no need to specify
# If counts.csv, then comma as delimiter. If counts.tsv, then \t as delimiter
# Ex: bash get_counts_columns.sh counts.tsv names.txt outputname.txt

OGCOUNTS=${1}
DELIM_NAME=${1##*.}
NAMES=${2}
OUTNAME=${3}

if [[ ${DELIM_NAME} == "csv" ]]; then
echo "Splitting using csv format"
DELIM=","
elif [[ ${DELIM_NAME} == "tsv" ]] || [[ ${DELIM_NAME} == "txt" ]]; then
echo "Splitting using tsv format"
DELIM="	"
else
echo "IDK what delimiter to use"
exit
fi

# First get counts header into column in a temp file
TMP=${OGCOUNTS}.tmp
for i in $(head -n 1 ${OGCOUNTS} | cut -f1- --output-delimiter "	"); do echo $i; done > ${TMP}

# Next find line number that matches query
#grep -n -f ${NAMES} ${TMP} | awk -F ":" '{print $1}' > ${TMP}1
echo "1" > ${TMP}1
for n in $(cat ${NAMES}); do grep -n ${n} ${TMP} | awk -F ":" '{print $1}'; done >> ${TMP}1

	# Sanity check. $NAMES should be a subset of the counts table headers
	# If I dont find all names in counts table, something is wrong
if [[ $(wc -l ${TMP}1 | awk '{print $1}') != $(( $(wc -l ${NAMES} | awk '{print $1}')+1 )) ]]; then
echo "${2} is not a subset of ${1}. Please check naming convention."
echo "$(wc -l ${TMP}1 | awk '{print $1}') / $(( $(wc -l ${NAMES} | awk '{print $1}')+1 ))"
#exit
fi

# Add 1 to line numbers bcuz lines have been shifted down by one
#awk '{print $1+1}' ${TMP}1 > ${TMP}2 && mv ${TMP}2 ${TMP}1

# Now print line numbers to string
#CUTLINES=$(for i in $(cat ${TMP}1); do printf "${i},"; done)
#CUTLINES=${CUTLINES%,}
#CUTLINE=$(printf "1,${CUTLINES}")
#echo $CUTLINE
counter=1
for i in $(cat ${TMP}1); do cut -f ${i} ${OGCOUNTS} > ${OGCOUNTS}.cut${counter}; counter=$(( $counter + 1 )); done

# Now get columns
#cut -f ${CUTLINE} --output-delimiter "	" ${OGCOUNTS} > ${OUTNAME}
\rm ${TMP} ${TMP}1
paste $(ls ${OGCOUNTS}.cut* | awk -F "cut" '{print $2, $0}' | \sort -n -k1 | sed 's/^[0-9][0-9]* //') > ${OUTNAME}
\rm ${OGCOUNTS}.cut*

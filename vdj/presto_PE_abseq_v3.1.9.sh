#!/usr/bin/env
#   Authors Alex Pan (apan82@gmail.com) and Altan Turkoglu
#   7/24/18

### Parses paired-end Abseq v3 fastqs through the pRESTO Abseq v3.1.9 Workflow
### Run in the directory containing fastq files
### Required inputs: 2 fastq files, R1 C primers fasta, R2 TS primers fasta,
### ..., internal C region fasta, and an immune reference fasta
### Also input # of cores to be used
### bash presto_PE_abseq_v3.1.9.sh R1.fastq R2.fastq ~osu9900/clonality_scripts/preMixcr/Abseq_scripts/AbSeqV3_Human_R1CPrimers.txt.fasta AbSeqV3_Human_R2TSPrimers.txt.fasta AbSeqV3_Human_InternalCRegion.txt.fasta 


# Load python 3.6
#module load python/3.6
module load python/3.6-conda5.2

# Assign arugments to variables
INPUTR1=${1}
INPUTR2=${2}
R1C=${3}
R2TS=${4}
INTC=${5}
IMREF=${6}
N_CORES=${7}

# Prefixes of outputs
TEMPR1=${INPUTR1%%_*}_R1
TEMPR2=${INPUTR2%%_*}_R2

#gunzip fastq files if necessary
if [ ${INPUTR1##*.} == "gz" ]; then
	echo "Unziping" ${INPUTR1}
	gunzip ${INPUTR1}
	INPUTR1=${INPUTR1%.gz}
fi

if [ ${INPUTR2##*.} == "gz" ]; then
	echo "Unzipping" ${INPUTR2}
	gunzip ${INPUTR2}
	INPUTR2=${INPUTR2%.gz}
fi

# Subsample if necessary
#SEED=$SECONDS
#seqtk sample -s $SEED ${INPUTR1} 500000 > tmp_${TEMPR1}_seqtk && mv \
#	tmp_${TEMPR1}_seqtk ${INPUTR1}
#seqtk sample -s $SEED ${INPUTR2} 500000 > tmp_${TEMPR2}_seqtk && mv \
#	tmp_${TEMPR2}_seqtk ${INPUTR2}

### Filterseq and MaskPrimers for each paired-end read
time FilterSeq.py quality -s ${INPUTR1} -q 20 --outname ${TEMPR1} --log FS1.log --nproc $N_CORES
time FilterSeq.py quality -s ${INPUTR2} -q 20 --outname ${TEMPR2} --log FS2.log --nproc $N_CORES
ParseLog.py -l FS1.log FS2.log -f ID QUALITY


time MaskPrimers.py score -s ${TEMPR1}_quality-pass.fastq -p ${R1C} \
	--start 0 --mode cut --maxerror 0.2 --outname ${TEMPR1} --log MP1.log --nproc $N_CORES
time MaskPrimers.py score -s ${TEMPR2}_quality-pass.fastq -p ${R2TS} \
	--start 17 --barcode --mode cut --maxerror 0.5 --outname ${TEMPR2} --log MP2.log --nproc $N_CORES
ParseLog.py -l MP1.log MP2.log -f ID PRIMER BARCODE ERROR 



### PairSeq
time PairSeq.py -1 ${TEMPR1}_primers-pass.fastq -2 ${TEMPR2}_primers-pass.fastq \
	--2f BARCODE --coord illumina

### BuildConsensus
time BuildConsensus.py -s ${TEMPR1}_primers-pass_pair-pass.fastq -n 1 -q 0 --freq 0.6 --bf BARCODE --pf PRIMER \
	--prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname ${TEMPR1} --log BC1.log --nproc $N_CORES
time BuildConsensus.py -s ${TEMPR2}_primers-pass_pair-pass.fastq -n 1 -q 0 --freq 0.6 --bf BARCODE --pf PRIMER \
	--prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname ${TEMPR2} --log BC2.log --nproc $N_CORES
ParseLog.py -l BC1.log BC2.log -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR

### PairSeq once more
time PairSeq.py -1 ${TEMPR1}_consensus-pass.fastq -2 ${TEMPR2}_consensus-pass.fastq \
	--coord presto

### Assesmble Pairs
time AssemblePairs.py sequential -1 ${TEMPR1}_consensus-pass_pair-pass.fastq \
	-2 ${TEMPR2}_consensus-pass_pair-pass.fastq -r ${IMREF} \
	--coord presto --rc tail --scanrev --1f CONSCOUNT --2f CONSCOUNT PRCONS \
	--alpha 1e-05 --maxerror 0.3 --minlen 8 --maxlen 1000 --aligner blastn \
	--minident 0.5 --evalue 1e-05 --maxhits 100 --outname ${TEMPR1%_R1}-C --log AP.log --nproc $N_CORES
ParseLog.py -l AP.log -f ID REFID LENGTH OVERLAP GAP ERROR IDENTITY

### FilterSeq again
time FilterSeq.py maskqual -s ${TEMPR1%_R1}-C_assemble-pass.fastq \
	-q 0 --outname ${TEMPR1%_R1}-C --log FS_mkql.log --nproc $N_CORES
ParseLog.py -l FS_mkql.log -f ID MASKED

### MaskPrimers
time MaskPrimers.py align -s ${TEMPR1%_R1}-C_maskqual-pass.fastq -p ${INTC} --maxlen 50 \
	--skiprc --gap 1 1 --mode cut --maxerror 0.2 --outname ${TEMPR1%_R1}-C --log MP3.log --nproc $N_CORES
ParseLog.py -l MP3.log -f ID PRIMER ERROR

### ParseHeaders
ParseHeaders.py rename -s ${TEMPR1%_R1}-C_primers-pass.fastq -f PRIMER \
	-k CREGION --act first --outname ${TEMPR1%_R1}-C_h1

### pRESTOR Abseq Report goes here*******

### ParseHeaders
ParseHeaders.py collapse -s ${TEMPR1%_R1}-C_h1_reheader.fastq \
	-f CONSCOUNT --act min --outname ${TEMPR1%_R1}-C_h2
ParseHeaders.py table -s ${TEMPR1%_R1}-C_h2_reheader.fastq \
	-f ID PRCONS CREGION CONSCOUNT --outname ${TEMPR1%_R1}-C_h3

### CollapseSeq
### Missing select by 1st seq ???
time CollapseSeq.py -s ${TEMPR1%_R1}-C_h2_reheader.fastq -n 0 --inner \
	--uf PRCONS CREGION --cf CONSCOUNT --act sum --keepmiss \
	--outname ${TEMPR1%_R1}-C 

#ParseHeaders
ParseHeaders.py table -s ${TEMPR1%_R1}-C_collapse-unique.fastq \
	-f ID PRCONS CREGION CONSCOUNT DUPCOUNT --outname ${TEMPR1%_R1}-C_h4

### Presto Parition goes here
time SplitSeq.py group -s ${TEMPR1%_R1}-C_collapse-unique.fastq -f CONSCOUNT --num 2 --outname ${TEMPR1%_R1}-C

###
ParseHeaders.py table -s ${TEMPR1%_R1}-C_atleast-2.fastq -f ID PRCONS CREGION CONSCOUNT DUPCOUNT

### gzip output
gzip ${TEMPR1%_R1}-C_atleast-2.fastq

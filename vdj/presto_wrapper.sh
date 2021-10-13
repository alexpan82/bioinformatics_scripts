for sample in $(cat samples.txt); do
echo '
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=10
#PBS -N presto_'${sample}'
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=002
#PBS -A PAS1374


date

cd $PBS_O_WORKDIR

pwd
### Make directory struc for what weldr expects
if [ ! -d output ]; then
        mkdir output
fi

if [ ! -d output/'${sample}' ]; then
        mkdir output/'${sample}'
fi

### Get FASTQ names
cd raw
R1FASTQ=$(ls '${sample}'*R1*fastq.gz)
R2FASTQ=$(ls '${sample}'*R2*fastq.gz)
cd ../output/'${sample}'
mkdir presto
cd presto
\rm raw_R1.fastq.gz raw_R2.fastq.gz
cp ../../../raw/${R1FASTQ} raw_R1.fastq.gz
cp ../../../raw/${R2FASTQ} raw_R2.fastq.gz

# PRESTO
REFLOC=~osu9900/clonality_scripts/preMixcr/Abseq_scripts
bash ../../../presto_PE_abseq_v3.1.9.sh raw_R1.fastq.gz raw_R2.fastq.gz \
${REFLOC}/AbSeqV3_Human_R1CPrimers.txt.fasta \
${REFLOC}/AbSeqV3_Human_R2TSPrimers.txt.fasta ${REFLOC}/AbSeqV3_Human_InternalCRegion.txt.fasta \
${REFLOC}/Immune_GRCh38.fasta 10

\rm raw_R1.fastq.gz raw_R2.fastq.gz

date
echo "Workflow complete"
' | qsub

done

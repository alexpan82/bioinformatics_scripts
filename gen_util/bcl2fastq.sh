#needs bcl directory as input
#Need to change the --use-bases-mask option based on which base pairs to read or ignore. That info can be found in the RunInfo.xml file

echo '
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=15
#PBS -l mem=12GB
#PBS -N bcl2fastq
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=002
#PBS -A PAS1374


date

cd $PBS_O_WORKDIR

pwd

source activate ~/anaconda2/envs/myenv2.7/

bcl2fastq -R ./'${1}' -o ./demultiplex_Pearlly --barcode-mismatches 0 --use-bases-mask Y151,I6N2,N8,Y151

date
echo "Workflow complete"
' | qsub

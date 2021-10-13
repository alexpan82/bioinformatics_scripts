echo '
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=10
#PBS -N posdup
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=002
#PBS -A PAS0472

date

cd $PBS_O_WORKDIR

pwd

source activate ~/anaconda2/envs/myenv2.7/

python ../bam_dup_v3.py -b '${1}' -o pos

date
echo "Workflow complete"
' | qsub

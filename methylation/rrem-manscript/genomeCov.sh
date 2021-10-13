echo '
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=10
#PBS -N job_wrapper
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=002
#PBS -A PAS0472


date

cd $PBS_O_WORKDIR
#source activate ~/anaconda2/envs/myenv2.7/
pwd

bedtools genomecov -bg -ibam '${1}' > '${1%bam}'genomeCov.bed

echo "Workflow complete"
' | qsub


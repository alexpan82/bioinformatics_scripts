echo '
#PBS -N bismark_genome_prep
#PBS -l nodes=1:ppn=18
#PBS -l walltime=4:00:00
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=022

cd $PBS_O_WORKDIR

pwd

/users/PAS0472/osu6369/tools/bismark_v0.17.0/bismark_genome_preparation --path_to_bowtie ~osu6369/tools/bowtie2-current/ --verbose '${1}'








' |qsub

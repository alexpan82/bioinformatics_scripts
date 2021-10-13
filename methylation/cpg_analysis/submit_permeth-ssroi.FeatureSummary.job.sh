#!/bin/bash

date=($( date '+%Y%m%d'))
echo '
#PBS -N job.'$date'_feat_summary
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00
#PBS -l mem=4GB
#PBS -S /bin/bash
#PBS -j oe

cd $PBS_O_WORKDIR
module load python/2.7.8

python /users/PAS0472/osu6369/jobs/bisulfite/analysis/YYMMDD_5x_permeth-ssroi.feature_summary.py '$date'

date
' | qsub 2>/dev/null


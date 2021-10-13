###### Input location of methylKit.R script you would like to run in "--------"

echo '#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=methylkit
#SBATCH --account=PASXXXX

date
pwd



cd $SLURM_SUBMIT_DIR
source activate ~osu9900/anaconda2/envs/r_env/
Rscript ./diffMeth_methylkit.R


date
echo "Workflow complete"
' | sbatch

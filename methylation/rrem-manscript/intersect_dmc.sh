# 1st argument is project to charge
# 2nd argument is DMCs in BED format
# 3rd argument is region of interests in BED format
echo '#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --job-name=intersect_'${3#*refgene.}'
#SBATCH --account='${1}'


date

cd $SLURM_SUBMIT_DIR

pwd

ls '${2}' '${3}' 

ROINAME='${3#*refgene.}'
ROINAME=${ROINAME%.hg38*}.bed

bedtools intersect -wo -bed -f 1E-9 -a '${2}' -b '${3}' > '${2%.*}'.${ROINAME}

date
echo "Workflow complete"
' | sbatch

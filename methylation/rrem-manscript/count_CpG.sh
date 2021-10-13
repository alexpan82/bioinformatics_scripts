echo '
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -N job_wrapper
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=002
#PBS -A PAS0472


date

cd $PBS_O_WORKDIR

pwd
#for i in '${1}'*bed; do echo -e "${i}\t$(cut -f 4 ${i} | \sort | uniq | wc -l)"; done
for i in '${1}'*bed; do echo -e "${i}\t$(cut -f 4,8,9 ${i} | \sort | uniq | awk -v n='${2}' '"'"'{if($2>=n){avg+=$3; total+=1}} END {print avg/total"\t"total}'"'"')"; done
date
echo "Workflow complete"
' | qsub

#!/bin/bash


###
##
#	MODIFIED TO WORK WITH PERMETH FILES
#	".permeth.gz" files should be in directory "../permeth/	
##
###


for samp in $(cat test.txt);do




#test
#for roi in 152*1kb.1kb*.roi;do
#all
for roi in *.roi;do



#
#check to see if .ssroi file exists
#

id=${samp##*/}
id=${id%.permeth}

dir=${roi%.roi}

if [ -s "$dir/$id"* ];then

	echo "skipping file "${id}"for feature "${dir}

else
	echo ${id}" ssroi "${dir}


echo '


#PBS -N job.permeth.ssroi_'$id'_'$dir'
#PBS -l nodes=1:ppn=2
#PBS -l mem=8gb
#PBS -l walltime=3:00:00
#PBS -S /bin/bash
#PBS -j oe
#PBS -W umask=022

cd $PBS_O_WORKDIR

date
SAMP='$samp'
ID='$id'
#ID=${SAMP%.sam*};
#ID=($( echo ${ID%.sam} ));
#ID=${ID##*/}
#
#make a short id to match the .pym file
#
#sID=($( echo ${ID%%.*}. ));


 
roi='$roi'
dir=($( echo ${roi%.roi} ));
if [ -d $dir ];then
	echo "output directory is "${dir}
else 
	mkdir ${dir}
	echo "Makeing output directory "${dir}
fi

echo "roi file is:"
ls ${roi}
echo "first line of "${roi}" is:"
head -2 ${roi} 
echo "Permeth file is:"
ls ${SAMP}

#load oakley python module
module load python/2.7.1

python /users/PAS0472/osu6369/jobs/bisulfite/analysis/ssroi.permeth.py ${SAMP} ${roi} ${ID}.${dir}.ssroi 5;
mv ${ID}.${dir}.ssroi ./${dir}/${ID}.${dir}.ssroi
echo "finished creating and moving "${ID}"."${dir}".ssroi"
#done

date

' | qsub 2>/dev/null

fi


done;
done;

#python ../../correlation_cov_IGNORE01.py -cutoff 100 -cpg1 '${1}' -cpg2 '${2}' > '${3}'
#for i in $(ls */*.txt | grep -v 'CpG.txt' | grep -v 'smooth'); do tail -n +4 ${i} | head -n 1000 > tmp && mv -f tmp ${i}; echo $i; done
#for i in $(ls */*.txt | grep -v 'CpG.txt' | grep -v 'smooth'); do n=${i%_*}; python ../smooth_expandingWindow.py -c 2 -t ${i} > ${n}_E_${i##*_}.smooth.txt; echo ${n}_E_${i##*_}.smooth.txt; done
#for i in $(ls */*.txt | grep -v 'CpG.txt' | grep -v 'smooth'); do n=${i%_*}; python ../smooth_expandingWindow.py -c 3 -t ${i} > ${n}_R_${i##*_}.smooth.txt; echo ${n}_R_${i##*_}.smooth.txt; done
#for i in */*smooth.txt; do awk '{print $1+1"\t"$2}' ${i} > tmp && mv -f tmp ${i}; done

#files=(correlation_cov_588N correlation_cov_608N_B4688 correlation_cov_608N_B4741 correlation_cov_617Na correlation_cov_617Nb correlation_cov_D101 correlation_cov_D4118)

#regions=(3putr 5putr cgi distalprom exons genebody intergenic introns proxprom rptmask)
#files=(correlation_cov_588N correlation_cov_D4118)
#regions=(cgi)
#regions=(epic)

for f in $(ls -d */); do
	cd ${f}
	efile1=$(ls *E_gauss.EPIC.txt.smooth.txt)
	efile2=$(ls *E_100.EPIC.txt.smooth.txt)
	rfile1=$(ls *R_gauss.EPIC.txt.smooth.txt)
	rfile2=$(ls *R_100.EPIC.txt.smooth.txt)
	list=($efile1 $efile2 $rfile1 $rfile2)
	for i in ${list[@]}; do
		#tail -n +3 ${i} | head -n 200 > ${i}.tmp
		head -n 200 ${i} > ${i}.tmp
	done
		
		#paste ${efile1}.tmp ${efile2}.tmp > tmpe
		#paste ${rfile1}.tmp ${rfile2}.tmp > tmpr
	paste ${efile1}.tmp ${efile2}.tmp ${rfile1}.tmp ${rfile2}.tmp > tmpER

	\rm *.tmp

	awk -F "\t" '{print $1,$2/$4,$5,$6/$8}' OFS="\t" tmpER > ${efile1%%_*}_ratio.EPIC.smooth.tsv
		#awk '{print $1"\t"$2/$4}' tmpr > ${rfile1%_gauss_*}_gauss_norm_100_ignore01.smooth.txt

		#\rm tmpe tmpr
	\rm tmpER

	echo ${efile1%%_*}_ratio.EPIC.smooth.tsv
	#done
	cd ..
done

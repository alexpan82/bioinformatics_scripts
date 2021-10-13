#!/bin/bash

### TRIMMING QC ###

###
#   Script to extract summary from all trimGalore trimming *trimming_report.txt files in located in trimmed reads subdirectories (*/*trimming_report.txt)
#   Matches for both paired-end and SR input sequences
#   David Frankhouser 2016-10-11
#
###


date=$(date +"%Y%m%d")



if [ -a "*.${date}tmp" ];then
    printf ":::ERROR:::\n\tThis script makes a number of \".${date}tmp\" files and then deletes them\n\tThis directory currently contains the following \".${date}tmp\" files:\n"
    ls *.${date}tmp
    printf "Remove these files before running!\n"
    exit 0
fi

outname=$1
if [ -z ${outname} ];then
    outname=${date}_BS-trimGaloreStats.txt
    printf ":::WARNING:::\n\tNo output name specified!\n\tUsing default output name: ${outname}\n"
    sleep 2s
fi

#make headers
read="R1"
printf "sampleName\n">name.${date}tmp
printf "adapterSeq(R1)\n">ad-${read}.${date}tmp
printf "ProcessedReads(R1)\n">pr-${read}.${date}tmp
printf "TrimmedReads(R1)\n">tr-${read}.${date}tmp
printf "50bpDetectedAdapter(R1)\n">ft-${read}.${date}tmp
printf "QualityTrimmed(R1)\n">qt-${read}.${date}tmp
printf "RRBS-2bpTrimmed(R1)\n">2bp-${read}.${date}tmp
read="R2"
printf "adapterSeq(R2)\n">ad-${read}.${date}tmp
printf "ProcessedReads(R2)\n">pr-${read}.${date}tmp
printf "TrimmedReads(R2)\n">tr-${read}.${date}tmp
printf "50bpDetectedAdapter(R2)\n">ft-${read}.${date}tmp
printf "QualityTrimmed(R2)\n">qt-${read}.${date}tmp
printf "RRBS-2bpTrimmed(R2)\n">2bp-${read}.${date}tmp
printf "TotalRemovedReads\n">tsr.${date}tmp
printf "RemovedReads(percent)\n">psr.${date}tmp

#Process files
printf "Processing bismark report files...\n"
for fname in fastq/*/*R1*trimming_report.txt;do
fbase=$(dirname ${fname})
printf "\t${fbase}:\n"
printf "${fbase}\n">>name.${date}tmp

for read in R1 R2;do
    file=$(ls ${fbase}/*${read}*trimming_report.txt) 
    printf "\t\t${read}...\n"
    grep $"Adapter sequence:" ${file} | awk '{print($3)}' >>ad-${read}.${date}tmp
    grep $"Processed reads:" ${file} | awk '{print($3)}' >>pr-${read}.${date}tmp
    grep $"Trimmed reads: " ${file} | awk '{print($3)}' >>tr-${read}.${date}tmp
    awk 'BEGIN{FS="\t"} {if ($1=="50") {print($2)} }' ${file}>>ft-${read}.${date}tmp
    grep $"Sequences were truncated to a varying degree because of deteriorating qualities" ${file} | awk 'BEGIN{FS="\t"} {print($2)}' >>qt-${read}.${date}tmp
    grep $"RRBS reads trimmed by 2 bp at the start when read started with CAA " ${file} | awk 'BEGIN{FS="\t"} {print($2)}' >>2bp-${read}.${date}tmp
    if [[ ${read} == "R2" ]];then
        grep $"Number of sequence pairs removed because at least one read was shorter than the length cutoff" ${file} | awk '{print($(NF-1))}' >>tsr.${date}tmp
        grep $"Number of sequence pairs removed because at least one read was shorter than the length cutoff" ${file} | awk '{print($NF)}' | sed 's/(//g' | sed 's/)//g' >>psr.${date}tmp
    fi
done



done

paste name.${date}tmp ad-R1.${date}tmp pr-R1.${date}tmp tr-R1.${date}tmp ft-R1.${date}tmp qt-R1.${date}tmp 2bp-R1.${date}tmp ad-R2.${date}tmp pr-R2.${date}tmp tr-R2.${date}tmp ft-R2.${date}tmp qt-R2.${date}tmp 2bp-R2.${date}tmp tsr.${date}tmp psr.${date}tmp>${outname}

printf "finished!\n"

#clean up
for i in *.${date}tmp;do rm ${i};done



### ALIGNMENT QC ###

###
#   Script to extract summary from all bismark report.txt files in located in directory
#   Matches for both paired-end and SR input sequences
#   David Frankhouser 2016-08-14
#
###


date=$(date +"%Y%m%d")



if [ -a "*.${date}tmp" ];then
    printf ":::ERROR:::\n\tThis script makes a number of \".${date}tmp\" files and then deletes them\n\tThis directory currently contains the following \".${date}tmp\" files:\n"
    ls *.${date}tmp
    printf "Remove these files before running!\n"
    exit 0
fi

outname=$1
if [ -z ${outname} ];then
    outname=${date}_BS-alignStats.txt
    printf ":::WARNING:::\n\tNo output name specified!\n\tUsing default output name: ${outname}"
    sleep 2s
fi

#make headers
printf "name\n">name.${date}tmp
printf "CHH-IC\n">chh.${date}tmp
printf "CHG-IC\n">chg.${date}tmp
printf "TotalReads\n">tot.${date}tmp
printf "UniqAlignedReads\n">aln.${date}tmp
printf "PercentAligned\n">per.${date}tmp
printf "unalignedFromMultiMapping\n">mm.${date}tmp
printf "expWideMeth\n">meth.${date}tmp

#process files
printf "Processing bismark report files...\n"
for file in bismark*/*report.txt;do
printf "\t${file}...\n"

printf "${file}\n" >>name.${date}tmp
grep $"C methylated in CHH context:" ${file} | awk 'BEGIN{FS="\t"} {print($2)}' >>chh.${date}tmp
grep $"C methylated in CHG context:" ${file} | awk 'BEGIN{FS="\t"} {print($2)}' >>chg.${date}tmp
grep $"Sequence.* analysed in total:" ${file} | awk 'BEGIN{FS="\t"} {print($2)}' >>tot.${date}tmp
grep $"Number of .* with a unique best hit.*:" ${file} | awk 'BEGIN{FS="\t"} {print($2)}' >>aln.${date}tmp
grep $"Mapping efficiency:" ${file} | awk 'BEGIN{FS="\t"} {print($2)}' >>per.${date}tmp
grep $"Sequence.* did not map uniquely:" ${file} | awk 'BEGIN{FS="\t"} {print($2)}' >>mm.${date}tmp
grep $"C methylated in CpG context:" ${file} | awk 'BEGIN{FS="\t"} {print($2)}' >>meth.${date}tmp

done

paste name.${date}tmp tot.${date}tmp aln.${date}tmp per.${date}tmp mm.${date}tmp chh.${date}tmp chg.${date}tmp meth.${date}tmp >${outname}

printf "finished!\n"

#clean up
for i in *.${date}tmp;do rm ${i};done

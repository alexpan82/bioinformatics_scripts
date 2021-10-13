
python ../../analysis_scripts/cpg_analysis/attach_genesym.py $1 ~/scratch/roi_files/mm10/refseq_roi/refgene.mm10.bed genebody test1.txt

python ../../analysis_scripts/cpg_analysis/attach_genesym.py test1.txt ~/scratch/roi_files/mm10/refseq_roi/refgene.10kb.1kb.prom.mm10.bed distal_prom test2.txt

python ../../analysis_scripts/cpg_analysis/attach_genesym.py test2.txt ~/scratch/roi_files/mm10/refseq_roi/refgene.1kb.1kb.prom.mm10.bed proximal_prom test3.txt

python ../../analysis_scripts/cpg_analysis/attach_genesym.py test3.txt ~/scratch/roi_files/mm10/refseq_roi/refgene.exons.mm10.bed exons test4.txt

python ../../analysis_scripts/cpg_analysis/attach_genesym.py test4.txt ~/scratch/roi_files/mm10/refseq_roi/refgene.introns.mm10.bed introns test5.txt

python ../../analysis_scripts/cpg_analysis/attach_genesym.py test5.txt ~/scratch/roi_files/mm10/refseq_roi/refgene.intergenic.exclude.10kb.prom.mm10.roi intergenic ${1%.txt}_genesAttached.txt

rm test*.txt

# Takes 2 CpG.txt files from methylKit and makes a union

from sys import argv
from collections import defaultdict

def readCpG(txtfile):
	returndic = defaultdict(lambda: ['', 0, 0])
	with open(txtfile, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			try:
				name = line[0]
				strand = line[3]
				coverage = int(line[4])
				numCs = round((float(line[5])/100.00)*coverage, 2)
				returndic[name] = [strand, coverage, numCs]
			except:
				pass
	return returndic

scripts, cpg1, cpg2 = argv

cpg1dic = readCpG(cpg1)
cpg2dic = readCpG(cpg2)

print('chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT')
for cpg in (set(cpg1dic.keys()) | set(cpg2dic.keys())):
	cov = cpg1dic[cpg][1] + cpg2dic[cpg][1]
	numCs = cpg1dic[cpg][2] + cpg2dic[cpg][2]
	strand = cpg1dic[cpg][0] + cpg2dic[cpg][0]
	tmp = cpg.split('.')
	chrm = tmp[0]
	pos = tmp[1]
	freqC = round(100*(numCs/cov), 2)
	freqT = round(100.0-freqC, 2)

	if len(strand) == 2:
		if strand[0] != strand[1]:
			continue
	print('\t'.join(list(map(lambda x: str(x), [cpg, chrm, pos, strand[0], cov, freqC, freqT]))))

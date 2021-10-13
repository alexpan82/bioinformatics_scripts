# Interleave only read pairs shared b/t PE fastq files

import gzip
from sys import argv
from collections import defaultdict

# Read fastq
def openfastq(r1):
	fastqdic = defaultdict(lambda:'')
	if r1.endswith('gz'):
		openedfile = gzip.open(r1, 'r')
	else:
		openedfile = open(r1, 'r')

	for n, line in enumerate(openedfile):
		if n % 4 == 0:
			name = line.strip().split()[0]
		fastqdic[name] += line
	return fastqdic

script, R1, R2 = argv

r1dic = openfastq(R1)
r2dic = openfastq(R2)

#outfile = open('interleave.fastq', 'w')

for name in r1dic.keys():
	if r2dic[name] != '':
		#outfile.write(r1dic[name]+r2dic[name])
		print(r1dic[name]+r2dic[name], end='')

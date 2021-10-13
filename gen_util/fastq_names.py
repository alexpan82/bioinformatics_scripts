# Grabs only fastq sequences with certain names

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
			name = line.strip().split()[0].lstrip('@')
		fastqdic[name] += line
	return fastqdic

R1 = argv[1]
R2 = argv[2]
namefile = argv[3]
try:
	outname = argv[4]
except:
	outname = 'subset'

r1dic = openfastq(R1)
r2dic = openfastq(R2)

nameset = set()
with open(namefile, 'r') as f:
	for line in f:
		line = line.strip()
		nameset.add(line)

outfiler1 = open(outname + '_R1.fastq', 'w')
outfiler2 = open(outname + '_R2.fastq', 'w')

for name in nameset:
	outfiler1.write(r1dic[name])
	outfiler2.write(r2dic[name])

outfiler1.close()
outfiler2.close()

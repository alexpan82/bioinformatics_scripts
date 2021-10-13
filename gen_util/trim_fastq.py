# Trims a set # of bp from the beginning or end of fastq
# Use '5' for trimming 5' end or '3' for 3'
# Ex: python trim_fastq.py R1.fastq.gz 17 5

import gzip
from sys import argv
from collections import defaultdict

# Read fastq
def openfastq(r1, trim, end):
	trim = int(trim)
	fastqdic = defaultdict(lambda:'')
	if r1.endswith('gz'):
		openedfile = gzip.open(r1, 'r')
	else:
		openedfile = open(r1, 'r')

	for n, line in enumerate(openedfile):
		line = line.strip()
		if (n+1) % 2 == 0:
			if end == '5': 
				print(line[trim:])
			elif end == '3':
				print(line[:-trim])
			else:
				quit()
		else:
			print(line)

script, R1, trim, end = argv

openfastq(R1, trim, end)

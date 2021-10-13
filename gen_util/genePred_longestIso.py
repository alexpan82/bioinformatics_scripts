# Takes a genePred formatted tsv
# Outpus only the longest isoform for each gene_id in genePred format to stdout

from collections import defaultdict
from sys import argv

script, inputFile = argv

# Hold length of each transcript
longest = defaultdict(lambda: 0)
longest_name = defaultdict(lambda: '')

f = open(inputFile, 'r')

for line in f:
	line = line.strip().split('\t')
	length = 0
	name1 = line[1]
	name2 = line[12]
	exon_starts = line[9].rstrip(',').split(',')
	exon_ends = line[10].rstrip(',').split(',')
	for start, end in zip(exon_starts, exon_ends):
		length += (int(end) - int(start))

	if length > longest[name2]:
		longest[name2] = length
		longest_name[name2] = name1
f.close()

f = open(inputFile, 'r')
for line in f:
	line = line.strip()
	l = line.split('\t')
	name1 = l[1]
	name2 = l[12]
	if longest_name[name2] == name1:
		print(line)
	
f.close()	



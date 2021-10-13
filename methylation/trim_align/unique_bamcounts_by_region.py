# Script to be used in conjuction with bam_overlaps_region.sh
# Bins unique bam reads by region name/class. For example, if a unique read is associated with rpt class x and y, the x and y bins get iterated by 1. If a unique read is associated with 2 elements that are both of class x, x only gets iterated by 1.
# This script is optimized for repeat elements so column changes are needed for other rois

from sys import argv
from collections import defaultdict

script, bed, out = argv

bed_line = open(bed, 'r')
outfile = open(out, 'w')
bed_line.readline()

# Create dictionary with this structure: {qname:{class}}
dic = defaultdict( lambda: set() )


for line in bed_line:
	line = line.strip().split('\t')
	try:
		name = line[3]
		genesym = line[8]
		dic[name].add(genesym)
	except:
		pass

bed_line.close()

# make dictionary with this structure: {class:{counts}}
histogram = defaultdict( lambda: 0 )

for qname in dic.keys():
	for i in dic[qname]:
		histogram[i] += 1

for n in sorted(histogram.keys()):
	outfile.write("%s\t\t\t%s\n" % (n, histogram[n]))

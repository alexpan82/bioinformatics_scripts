# Makes histogram of the number of times reads in an ambig_bam file occurs

import pysam
from sys import argv
from collections import defaultdict
import matplotlib.pyplot as plt

script, bam = argv

# Create dictionary to hold names and # of occurances 

names = defaultdict(lambda: 0)

with pysam.AlignmentFile(bam, "rb") as f:
	for line in f:
		names[line.qname] += 1

hist = []

for i in names.keys():
	hist.append(names[i])

# Plot histogram
binlist = range(1, len(set(hist))+1)
plt.hist(hist, bins=binlist)
plt.title("# of multimapping locations per unique read")
plt.xlabel("# of multimapping locations")
plt.ylabel("# of Unique reads")

fig = plt.show()

### This script will make a histogram of # of unique multi-mapped reads vs the # of different 
### families the read multimapped to
### For Ex: an amib read multimapped to AluY and AluSx. It gets binned in column "2"
### Takes a bed file from bedtools intersect with read name and annotation in colmn 4 & 9
### Histogram to also displays stacked bar for a annotation(s) of a 
#specific pattern (separated into all one annotation, none, and mixed)
	### Ex: To display the same info but only for reads that only get annotated as Alu*
	# set -a Alu


### Written by: Alex Pan (apan82@gmail.com)
import argparse
import re
from collections import defaultdict
from sys import argv
import csv
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

# If every element of a list starts with the pattern, returns 1. If some elements start with patt, return 0. If none, -1
def strsearch(lst, pat):
	count = 0
	for i in lst:
		if i.startswith(pat):
			count += 1
	if count == len(lst):
		return int(1)
	elif count < len(lst) and count > 0:
		return int(0)
	else:
		return int(-1)


def main(args):

	# Parse arguments
	bed = args.bed
	patt = args.patt

	# Declare dictionary to hold values
	dic = defaultdict(lambda: [])

	with open(bed, 'r') as f:	
		reader = csv.reader(f, delimiter='\t')	
		for line in reader:
			dic[line[3]].append(line[8])
		f.close()

	# Make list to hold histogram values
	totalhist = []
	pattall = []
	pattmix = []
	pattnone = []

	# If user specified pattern...

	for name in dic.keys():
		tmpset = set(dic[name])
		totalhist.append(len(tmpset))		
		if strsearch(tmpset, patt) == 1:
			pattall.append(len(tmpset))
		elif strsearch(tmpset, patt) == 0:
			pattmix.append(len(tmpset))
		elif strsearch(tmpset, patt) == -1:
			pattnone.append(len(tmpset))
	f.close()




	# Make stacked histogram. Assumes bed file is the complete set and annotations
	# that match the pattern and/or are in the optional file are a SUBSET of the complete set


	tmppattall = np.histogram( pattall, bins=max(totalhist)+1, range=(1,max(totalhist)+1) )
	tmppattmix = np.histogram( pattmix, bins=max(totalhist)+1, range=(1,max(totalhist)+1) )
	tmppattnone = np.histogram( pattnone, bins=max(totalhist)+1, range=(1,max(totalhist)+1) )
	data = (tmppattall[0], tmppattmix[0], tmppattnone[0],) 
	legend = ('All'+patt, 'Mixed'+patt, 'None'+patt,)

	print data[0]
	print data[1]
	print data[2]
	# Plot histogram
	binlist = range(1, len(data[0])+1)


	plt.bar(binlist, data[0], color='g')
	plt.bar(binlist, data[1], bottom=data[0], color='r')
	plt.bar(binlist, data[2], bottom=[x + y for x, y in zip(data[0], data[1])])

	
	plt.legend(legend)
	plt.xlabel("# of Multimapped Annotations ")
	plt.ylabel("# of Unique Multimapped Reads")

	fig = plt.show()
	



if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Removes primers from paired-end Ab-seq fastqs')
        parser.add_argument('-b', dest='bed', type=str, help='bed file from bedtools intersect with read name and family annotation in colmns 4 & 9')
        parser.add_argument('-a', dest='patt', default = "no", type=str, help='Optional: Type pattern (If you want the histogram to also display separate bars for a annotation(s) of a specific pattern)')
	args=parser.parse_args()
        main(args)


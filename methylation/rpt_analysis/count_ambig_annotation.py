### This script will make a histogram of # of unique multi-mapped reads vs the # of different 
### families the read multimapped to
### For Ex: an amib read multimapped to AluY and AluSx. It gets binned in column "2"
### Takes a bed file from bedtools intersect with read name and annotation in colmn 4 & 9
### Optional arguements:
	### If you want the histogram to also display separate bars for a annotation(s) of a 
	#specific pattern
		### Ex: To display the same info but only for reads that only get annotated as Alu*
		# set -a Alu
	### If you also want to display separate bars for specific reads, please pass a file
	# containing read names (1 per line) to -f

### Written by: Alex Pan (apan82@gmail.com)

import argparse
import re
from collections import defaultdict
from sys import argv
import csv
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

# If every element of a list starts with the pattern, returns true
def strsearch(lst, pat):
	tmp = 1
	for i in lst:
		if i.startswith(pat):
			tmp = 1
		else:
			tmp = 0
			break
	return tmp


def main(args):

	# Parse arguments
	bed = args.bed
	patt = args.patt
	xfile = args.xfile

	# Declare dictionary to hold values
	dic = defaultdict(lambda: [])

	with open(bed, 'r') as f:	
		reader = csv.reader(f, delimiter='\t')	
		for line in reader:
			dic[line[3]].append(line[8])
		f.close()

	# Make list to hold histogram values
	totalhist = []
	patthist = []

	# If user specified pattern...
	if patt != "no":
		for name in dic.keys():
		#tmp = float(stats.mode(dic[name])[1][0])/float(len(dic[name]))
			tmpset = set(dic[name])
			totalhist.append(len(tmpset))		
			if strsearch(tmpset, patt) == 1:
				patthist.append(len(tmpset))
			else:
				pass
		f.close()
	else:
		for name in dic.keys():
			totalhist.append(len(set(dic[name])))


	# If user specified read name file, print to separate list
	xhist = []
	if xfile != "no":
		with open(xfile, 'r') as f:
			for name in f:
				name = name.strip()
				xhist.append(len(set(dic[name])))

			f.close()

	

	# Make stacked histogram. Assumes bed file is the complete set and annotations
	# that match the pattern and/or are in the optional file are a SUBSET of the complete set

	data = (totalhist,)
	legend = (bed,)
	if patt != "no":
		del data
		tmphist = np.histogram( totalhist, bins=max(totalhist)+1, range=(1,max(totalhist)+1) )
		tmppatt = np.histogram( patthist, bins=max(totalhist)+1, range=(1,max(totalhist)+1) )
		subhist = tmphist[0] - tmppatt[0]
		data = (subhist, tmppatt[0],) 
		legend = legend + (patt,)

	if xfile != "no":
		del data
		tmpxhist = np.histogram( xhist, bins=max(totalhist)+1, range=(1,max(totalhist)+1) )
		if patt != "no":
			# Here we assume xfile is a subset of the pattern histogram tmppatt
			subpatt = tmppatt[0] - tmpxhist[0]
			data = (subhist, subpatt, tmpxhist[0],)
		else:
			tmphist = np.histogram( totalhist, bins=max(totalhist)+1, range=(1,max(totalhist)+1) )
			subhist = tmphist[0] - tmpxhist[0]
			data = (subhist, tmpxhist[0],)
		legend = legend + (xfile,)

	
	# Plot histogram
	binlist = range(1, len(data[0])+1)

	if patt != "no" and xfile != "no":
		plt.bar(binlist, data[1], color='g')
		plt.bar(binlist, data[2], bottom=data[1], color='r')
		plt.bar(binlist, data[0], bottom=[x + y for x, y in zip(data[1], data[2])])

	elif patt == "no" and xfile != "no":
		plt.bar(binlist, data[2], color='r')
		plt.bar(binlist, data[0], bottom=data[2])

	elif patt != "no" and xfile == "no":
		plt.bar(binlist, data[1], color='g')
		plt.bar(binlist, data[0], bottom=data[1])

	else:
		plt.bar(binlist, data[0])

	
	plt.legend(legend)
	plt.xlabel("# of Multimapped Annotations ")
	plt.ylabel("# of Unique Multimapped Reads")

	fig = plt.show()
	



if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Removes primers from paired-end Ab-seq fastqs')
        parser.add_argument('-b', dest='bed', type=str, help='bed file from bedtools intersect with read name and family annotation in colmns 4 & 9')
        parser.add_argument('-a', dest='patt', default = "no", type=str, help='Optional: Type pattern (If you want the histogram to also display separate bars for a annotation(s) of a specific pattern)')
        parser.add_argument('-f', dest='xfile', default = "no", type=str, help='Optional: If you also want to display separate bars for specific reads, please pass a file containing read names (1 per line)')
	args=parser.parse_args()
        main(args)


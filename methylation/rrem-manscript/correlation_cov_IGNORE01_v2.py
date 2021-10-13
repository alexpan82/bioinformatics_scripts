# This script takes 2 CpG.txt files (methylKit)
# Will calculate correlation and output histograms of coverage for shared correlated CpGs for each file.
# Will then output the ratio of correlated Cpgs b/t each file at each coverage

# Ex: python correlation_coy.py -cutoff 10 -cpg1 [FILE1] -cpg2 [FILE2]

import argparse
import pandas as pd
from collections import defaultdict
import math
#import numpy as np
import scipy
from scipy import stats

# Open and read CpG.txt file and return dic
# dic[CpG] = [cov, meth]
def readCpG(txtfile):
	returndic = defaultdict(lambda: [0,0])
	with open(permeth1, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			try:
				name = line[0]
				coverage = int(line[4])
				meth = float(line[5])
				returndic[name] = [coverage, meth]
			except:
				pass
	return returndic

def main(args):

	# Read args
	if args.cutoff == 'e':
		cutoff = args.cutoff
	else:
		cutoff = float(args.cutoff)

	methobject1 = readCpG(args.cpg1)
	methobject2 = readCpG(args.cpg2)


	# Union of CpGs
	allCpGs = methobject1.keys() & methobject2.keys()

	# pair all CpG info b/t 2 files
	# methcov[name] = [[cov1, meth1],...]
	methcov = defaultdict(lambda: [])

	for c in allCpGs:
		if methobject1[c] != [0,0] and methobject2[c] != [0,0]:
			methcov[0].append(methobject1[c])
			methcov[1].append(methobject2[c])



	# Now plot coverage histogram for every sample vs every other sample
		# I put 2 strings in a set every loop so that I do not plot the same comparison twice
		# Prob a more elegent way to do this, but eh

	i = 0
	checkset = set()
	for sample2 in methcov.keys():
		checkstring1 = '%s, %s' % (sample1, sample2)

		abins = [0]*(maxcov+1)
		# Calculate pearson coefficient
		sample1meth = []
		sample2meth = []
		for cpg1,cpg2 in zip(methcov[sample1], methcov[sample2]):
			if sample1 in samplegroup[groupid]:
				cov1 = cpg1[0]
				meth1 = cpg1[1]
				meth2 = cpg2[1]

			else:
				cov1 = cpg2[0]
				meth1 = cpg2[1]
				meth2 = cpg1[1]
				
			# if methylation difference is less than cutoff, plot the cov of the sample with the 
			# user specified groupid
			if cutoff == 'e':
				# TEMPORARILY applying additional filter here
				if (meth1 == 1 and meth2 == 1) or (meth1 == 0 and meth2 == 0):
					pass
				else:
					abins[cov1] += math.exp(-0.1*abs(meth1-meth2)*abs(meth1-meth2)*100*100)
					# appened methlyation to calc pearson correlation
					sample1meth.append(cpg1[1])
					sample2meth.append(cpg2[1])
					
			else:
				if abs(meth1-meth2) <= cutoff:
					# TEMPORARILY applying additional filter here
					if (meth1 == 1 and meth2 == 1) or (meth1 == 0 and meth2 == 0):
						pass
					else:
						abins[cov1] += 1
						# appened methlyation to calc pearson correlation
						sample1meth.append(cpg1[1])
						sample2meth.append(cpg2[1])
				else:
					pass

			i += 1
			#print(abins)
			print(scipy.stats.pearsonr(sample1meth, sample2meth))
			if checkstring1 == 'Emseq, RRBS' or checkstring1 == 'RRBS, Emseq':
				for n,i in enumerate(abins):
					print("%s\t%s" % (n,i))
		checkstring2 = '%s, %s' % (sample2, sample1)
		checkset.add(checkstring1)
		checkset.add(checkstring2)



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='See comments in code')
	parser.add_argument('-cpg1', dest='cpg1', type=str, help='CpG.txt from methlKit')
	parser.add_argument('-cpg2', dest='cpg2', type=str, help='CpG.txt from methlKit')
	parser.add_argument('-cutoff', dest='cutoff', type=str, help='methylation difference cutoff (%). If -cutoff e, then will use gaussian binning and there will not be a discrete cutoff')
	args=parser.parse_args()
	main(args)

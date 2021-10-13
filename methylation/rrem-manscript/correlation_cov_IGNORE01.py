# This script takes a meth unite object from methlKit
# will output coverage distributions for every sample comparison whose CpGs meet a give CpG methylation difference cutoff
# Also give a tsv that lists sample names on each new line along with a group id in another column
	# Make sure to give names in the same order as what was used to generate the methlKit object
	# Make sure to give unique group ids in the file
# Ex: python correlation_coy.py -cutoff 10 -names names.txt -meth meth.tsv -id groupname

import argparse
import pandas as pd
from collections import defaultdict
import math
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
#import matplotlib.backends.backend_pdf
#import numpy as np
import scipy
from scipy import stats

def main(args):

	# Read args
	if args.cutoff == 'e':
		cutoff = args.cutoff
	else:
		cutoff = float(args.cutoff)/100.0

	meth = args.meth
	names = args.names
	groupid = args.groupid

	### Write to SAME pdf
	pngname = meth.rstrip('.tsv') + "_" + names.rstrip('.txt') + "." + str(cutoff) + "." + groupid +".png"

	# Read names and keep them in order
	# nameorder[order] = name
	# samplegroup[group] = [samplenames]
	nameorder = defaultdict(lambda: 'na')
	samplegroup = defaultdict(lambda: [])
	with open(names, 'r') as f:
		for n,line in enumerate(f):
			line = line.strip().split('\t')
			nameorder[n+1] = line[0]
			samplegroup[line[1]].append(line[0])


	# Read thru meth and associate data with sample names
	# methcov[name] = [[cov1, meth1],...]
	methcov = defaultdict(lambda: [])
	methoject = pd.read_csv(meth, sep='\t')

	for order in nameorder.keys():
		covname = 'coverage' + str(order)
		methname = 'numCs' + str(order)
		for c,m in zip(methoject[covname].tolist(), methoject[methname].tolist()):
			tmplist = [int(c), float(m)/float(c)]
			methcov[nameorder[order]].append(tmplist)


	# Decide how many subplots to make on one pdf page
	square = map(lambda x: x**2, range(1,21))
	# Find total number of comparisons b/t groups
	totalcomp = 0
	netcomp = 0
	for g in samplegroup.keys():
		for s in samplegroup[g]:
			totalcomp += 1
	for g in samplegroup.keys():
		tmp = 0
		for s in samplegroup[g]:
			tmp += 1
		totalcomp -= tmp
		netcomp += tmp * totalcomp

	for num in sorted(square):
		if num >= netcomp:
		        size = math.sqrt(num)
		        break
	# Find largest cov to bound histogram plots
	covset = set()
	for n in methcov.keys():
		for cpg in methcov[n]:
			covset.add(cpg[0])
	
	maxcov = max(covset)
	del covset
	# Now plot coverage histogram for every sample vs every other sample
		# I put 2 strings in a set every loop so that I do not plot the same comparison twice
		# Prob a more elegent way to do this, but eh

	#plt.rcParams.update({'font.size': 4})
	#mpl.rc('ytick', labelsize=4)
	#mpl.rc('xtick', labelsize=4)
        # Create subplots
	#fig = plt.figure()
	#plt.axis('off')
	#ax = {}
	#for n in range(0, netcomp):
	#	ax[n+1] = fig.add_subplot(size, size, n+1)
	i = 0
	checkset = set()
	for sample1 in methcov.keys():
		for sample2 in methcov.keys():
			checkstring1 = '%s, %s' % (sample1, sample2)
			# plot empty plot if comparing the same sample
			if sample1 == sample2:
				pass
			# empty plot if comparison has already been done
			elif checkstring1 in checkset:
				pass
			# empty plot if comparison is within same group
			elif sample1 in samplegroup[groupid] and sample2 in samplegroup[groupid]:
				pass
			elif sample1 not in samplegroup[groupid] and sample2 not in samplegroup[groupid]:
				pass
			else:
				#if sample1[0:4] in sample2:
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
				'''
				ax[i+1].bar(np.arange(len(abins)), abins)
				ax[i+1].xaxis.set_major_locator(plt.MaxNLocator(4))
				ax[i+1].tick_params(axis='y', direction='in', length=1)
				ax[i+1].tick_params(axis='x', length=1)
				ax[i+1].yaxis.set_major_locator(plt.MaxNLocator(5))
				ax[i+1].set_xlim(left=0, right=maxcov)
				ax[i+1].set_title(checkstring1)
				'''
				i += 1
				#print(abins)
				print(scipy.stats.pearsonr(sample1meth, sample2meth))
				if checkstring1 == 'Emseq, RRBS' or checkstring1 == 'RRBS, Emseq':
					for n,i in enumerate(abins):
						print("%s\t%s" % (n,i))
			checkstring2 = '%s, %s' % (sample2, sample1)
			checkset.add(checkstring1)
			checkset.add(checkstring2)

	#print(i)
	#print(netcomp)
	#plt.tight_layout()
	#plt.savefig(pngname, format='png')
	#plt.close()



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='See comments in code')
	parser.add_argument('-meth', dest='meth', type=str, help='meth.tsv from methlKit unite object')
	parser.add_argument('-cutoff', dest='cutoff', type=str, help='methylation difference cutoff (%). If -cutoff e, then will uses gaussian binning and there will not be a hard cutoff')
	parser.add_argument('-names', dest='names', type=str, help='make sure to give names in the same order as what was used to generate the methlKit object')
	parser.add_argument('-id', dest='groupid', type=str, help='group id to plot coverage histogram for. must be one of the group ids specified in the file passed to -names')
	args=parser.parse_args()
	main(args)

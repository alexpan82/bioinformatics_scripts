# This script plots clustered barplot for the normalized number of Vgene usage across samples
	# Normalized by # of reads that aligned to vgenes per sample

# Input file is a design tsv
	# 1st column are file names (or paths to files)
	# 2nd column are file names

# Files referrenced by design.tsv are from mixcr exportAlignmixcr exportAlignments -f --preset full 
	#-nMutationsRelative V5UTR VTranscript \
	#-nMutationsRelative L1 VTranscript \
	#-nMutationsRelative L2 VTranscript \
	#-nMutationsRelative FR1 VTranscript \
	#-nMutationsRelative CDR1 VTranscript \
	#-nMutationsRelative FR2 VTranscript \
	#-nMutationsRelative CDR2 VTranscript \
	#-nMutationsRelative FR3 VTranscript \
	#-nMutationsRelative VCDR3Part VTranscript \
	#alignments.vdjca alignments.txt
# Script is written for mixcr v3.0.5 output. May need to be adjusted for other mixcr versions

# Counts the gene that had the highest alignment score
	# If more than one gene had the same score, those genes will each be counted once

# Written by: Alex Pan

import math
import numpy as np
from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import argparse

def Average(lst): 
    return sum(lst) / len(lst) 

def main(args):

	# parse args
	design = args.design
	celltype = args.celltype
	
	if celltype not in ['TRAV', 'TRBV', 'TRDV', 'TRGV', 'IGHV', 'IGKV', 'IGLV']:
		print 'Type TR[A,B,D,G]V or IG[H,K,L]V'
		quit()

	### Write to SAME pdf
	outname = '%sUsage.mixcr.pdf' % celltype.lower()
	pdf = PdfPages(outname)

	# Read from design file and store into array
	designdic = []
	with open(design, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			fname1 = line[0].rstrip('.align.txt')
			fname2 = line[1].rstrip('.align.txt')
			designdic.append([fname1, fname2])

	# vdic[sample name][vgene] = # of times vgene was used per sample
	vdic = defaultdict(lambda: defaultdict(lambda: 0))
	# Hold the number of reads that cover the vgene
	# norm[samplename][vgene] = coverage
	norm = defaultdict(lambda: 0)
	
	print 'Reading file...'
	### Read alignment files data table and store into a dictionary
	for pair in designdic:
		for fpath in pair:
			fpath = fpath + '.align.txt'
			with open(fpath, 'r') as f:
				for line in f:
					line = line.strip().split('\t')
					try:
						sample = fpath.rstrip('.align.txt')
						vgenelist = []
						qualitylist = []
						# Find the vgenes that a sequence aligned to
						# line[2] should correspond with 'allVHitsWithScore' column
						# For now, I will only consider the first gene in line[2]
							# It has the highest alignment score
						gname = line[2].split(',')[0]
						gname = gname.split('*')[0]
						if gname.startswith(celltype):
							vdic[sample][gname] += 1
							norm[sample] += 1
					except:
						pass

	print 'Normalizing... '
	# Lets normalize (# mutations / total read length for gene)
	normvdic = defaultdict(lambda: defaultdict(lambda: 0))
	for name in vdic.keys():
		for g in vdic[name].keys():
			normvalue = (float(vdic[name][g]) / float(norm[name]))
			normvdic[name][g] = normvalue

	del vdic


	# Decide how many subplots to make on one pdf page
	square = map(lambda x: x**2, range(1,21))
	for num in sorted(square):
		if num >= len(normvdic.keys()) / 2:
			size = int(math.sqrt(num))
			break

	### plot only the 1st ten vgenes. Throw everything else into 'Others'
	for f in normvdic.keys():
		other = 0.0
		popgenes = []
		#print f
		for index, gene in enumerate(sorted(normvdic[f].keys(), key=lambda u: normvdic[f][u], reverse=True)):
			if index > 10:
				other += normvdic[f][gene]
				popgenes.append(gene)
		if other != 0:
			normvdic[f]['Others'] = other
		for p in popgenes:
			normvdic[f].pop(p)
			
	#print normvdic
		# Create subplots
	print 'Creating subplots... '
	# Create subplots
	ax = {}
	fig = plt.figure()


	plt.rcParams.update({'font.size': 4})
	mpl.rc('ytick', labelsize=4)
	mpl.rc('xtick', labelsize=4)
	plt.axis('off')

	for n in range(0, len(normvdic.keys()) / 2):
		ax[n+1] = fig.add_subplot(size, size, n+1)

	# plot subplots
	for n, array in enumerate(designdic):
		# Set tickmarks for paired files
		xtickmarks = []
		for f in array:
			xtickmarks.append(f)

		# Assign colors to bars
		genecolors = {}
		for f in array:
			for g in normvdic[f].keys():
				genecolors[g] = np.random.rand(1,3)[0]
		# Make 'Others' gray
			genecolors['Others'] = np.array([0.1843, 0.3098, 0.3098])

		# Plot stacked bars
		# Also make legend
		# legend[gene] = (bar[x][0], frequncy)
		legendmark = {}

		for k, f in enumerate(array):
			bar = {}
			if 'Others' in normvdic[f].keys():
				data = [0] * len(array)
				data[k] = normvdic[f]['Others']
				# Color 'Others' bar first
				bar[0] = ax[n+1].bar(np.arange(len(data)), data, width=0.5, linewidth=0.1, color=genecolors['Others'])
				legendmark['Others'] = (bar[0][0], normvdic[f]['Others'])
				nextbot = normvdic[f]['Others']
			else:
				nextbot = 0
			# Color the rest of the stacked bars
			for index, g in enumerate(sorted(normvdic[f].keys(), key=lambda u: normvdic[f][u])):
				data = [0] * len(array)
				data[k] = normvdic[f][g]
				if g != 'Others':
					print np.arange(len(array))
					print data
					print genecolors[g]
					print nextbot
					bar[index+1] = ax[n+1].bar(np.arange(len(array)), data, width=0.5, linewidth=0.1, color=genecolors[g], bottom=nextbot)
					legendmark[g] = (bar[index+1][0], normvdic[f][g])
					nextbot += normvdic[f][g]


		# set other bar parameters
		ax[n+1].set_xticks(np.arange(len(array)))
		xtickNames = ax[n+1].set_xticklabels(xtickmarks)
		plt.setp(xtickNames, rotation=0, fontsize=4)
		ax[n+1].tick_params(axis='y', direction='in', length=1)
		ax[n+1].tick_params(axis='x', length=0)
		ax[n+1].yaxis.set_major_locator(plt.MaxNLocator(5))
		ax[n+1].yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%s%%' %(int(x*100.0))))
		#ax[n+1].set_title('%s (%s%%)     %s (%s%%)' % (name1, mut1, name2, mut2))
		ax[n+1].set_ylim(bottom=0, top=1)
		#ax[n+1].axhline(y=0.02)
		# Make legend
		legendtuple = tuple()
		nametuple = tuple()
		for g in sorted(legendmark.keys(), key=lambda u: legendmark[u][1], reverse=True):
			if g != 'Others':
				legendtuple += (legendmark[g][0],)
				nametuple += (g,)
		if 'Others' in normvdic[f].keys():
			legendtuple += (legendmark['Others'][0],)
			nametuple += ('Others',)
			ax[n+1].legend(legendtuple, nametuple, loc=0, fontsize=1.5)




	plt.tight_layout()
	plt.savefig(pdf, format='pdf')
	plt.close()


	pdf.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='See comments in code')
	parser.add_argument('-d', dest='design', type=str, help='Takes design tsv (See comments in code)')
	parser.add_argument('-c', dest='celltype', type=str, help='Type TR[A,B,D,G]V or IG[H,K,L]V')
	args=parser.parse_args()
	main(args)


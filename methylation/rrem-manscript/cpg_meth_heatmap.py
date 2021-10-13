# Script takes two CpG files (outputed from methylKit)
# Different in that it uses matplotlib to make heatmap
	# see https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html 
# Graphs the methylation of shared CpGs in a heatmap

import argparse
from collections import defaultdict
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def main(args):

	permeth1 = args.permeth1
	permeth2 = args.permeth2
	covdiff = float(args.covdiff)
	print(covdiff)
	pdfname = 'output_cov_%s.pdf' % int(covdiff)
	pdf = PdfPages(pdfname)

	# Create dictionary to store coverage info. Also insures we keep the data in order
	# dic['chromosome \t cpg position'] = coverage
	covdic1 = defaultdict(lambda: 0)
	covdic2 = defaultdict(lambda: 0)
	# Also create dictionaries to store methylation
	# dic['chromosome \t cpg position'] = methylation
	methdic1 = defaultdict(lambda: 'NA')
	methdic2 = defaultdict(lambda: 'NA')

	# Open and read
	with open(permeth1, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			try:
				name = line[0]
				coverage = int(line[4])
				meth = round(float(line[5]), 0)
				covdic1[name] = coverage
				methdic1[name] = meth
			except:
				pass
		f.close()

	with open(permeth2, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			try:
				name = line[0]
				coverage = int(line[4])
				meth = round(float(line[5]), 0)
				covdic2[name] = coverage
				methdic2[name] = meth
			except:
				pass
		f.close()


	# Put data into a numpy array 
	arrayoflists = []

		# Populate arrayoflists with [0,0,0,...]
		# Number of these arrays depend on unique cov of y (permeth2)
		# Number of 0's depends on unique cov of x (permeth1)

	for i in range(0, 101):
		# This is done because we use log scale to view
		tmp = [1]*(101)
		arrayoflists.append(tmp)
	arrayoflists = np.array(arrayoflists)

	# Populate arrayoflists
	# 1st make set with all CpGs between both files
	totalcpgs = set(covdic1.keys()).union(set(covdic2.keys()))
	#print(len(totalcpgs))
	countshared = 0
	for pos in totalcpgs:
		cov1 = covdic1[pos]
		cov2 = covdic2[pos]
		meth1 = methdic1[pos]
		meth2 = methdic2[pos]
		
		# Conditionals for coverage difference
		if meth1 == 'NA' or meth2 == 'NA':
			continue
		elif cov1 < covdiff and cov2 < covdiff:
			continue
		else:
			pass

		y = int(meth2)
		x = int(meth1)
		arrayoflists[100-y][x] += 1
		countshared += 1
	print(countshared)

	# Also keep track of max values
	maxvar = 0
	for y in arrayoflists:
		for x in y:
			if x > maxvar:
				maxvar = x
	#print arrayoflists
	matplotlib.rc('ytick', labelsize=6)
	matplotlib.rc('xtick', labelsize=6)

	# Make labels:
	xlabels = []
	ylabels = []
	for i in range(0,101):
		xlabels.append(i)
	for i in range(0, 101):
		ylabels.append(i)
	xlabels = sorted(xlabels)
	xlabels = list(map(str, xlabels))	
	ylabels = sorted(ylabels, reverse=True)
	ylabels = list(map(str, ylabels))
	# Plot heatmap
	fig, ax = plt.subplots()
	#im = ax.imshow(arrayoflists)
	# Log norm
	im = ax.pcolor(arrayoflists[::-1], norm=LogNorm(vmin=1, vmax=maxvar), cmap = 'turbo')
	cbar = ax.figure.colorbar(im, ax=ax, **{})
	cbar.ax.set_ylabel('', rotation=-90, va="bottom")
	ax.set_xticks(list(map(lambda x: x+0.5, np.arange(len(xlabels)))))
	ax.set_yticks(list(map(lambda y: y+0.5, np.arange(len(ylabels)))))
	#ax.set_ylim(top=ybins-min(set(covdic2.values()))+1)
	ax.set_xticklabels(xlabels)
	ax.set_yticklabels(ylabels[::-1])
	# Rotate labels
	plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
	# # in each box
	'''for i in range(len(ylabels)):
	    for j in range(len(xlabels)):
		text = ax.text(j, i, arrayoflists[i, j],
		               ha="center", va="center", color="w")
	'''
	ax.set_title("Saved CpGs")
	fig.tight_layout()
	plt.savefig(pdf, format='pdf')
	plt.close()
	pdf.close()
	
	

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='See comments in code')
        parser.add_argument('-x', dest='permeth1', type=str, help='Permeth file whose coverage to graph on the x axis')
        parser.add_argument('-y', dest='permeth2', type=str, help='Permeth file whose coverage to graph on the y axis')
        parser.add_argument('-cov', dest='covdiff', type=str, default='0', help='Coverage cutoff. For example: -diff 5 will only plot shared CpGs that have a coverage greater than or equal to 5')
        args=parser.parse_args()
        main(args)





# Script takes two CpG files (outputed from methylKit)
# Different in that it uses matplotlib to make heatplot
	# see https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html 
# Graphs the coverage of both in a heatplot

import argparse
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def main(args):

	pdf = PdfPages('output_meth50cov.pdf')
	permeth1 = args.permeth1
	permeth2 = args.permeth2
	ybins = 100

	# Create dictionary to store coverage info. Also insures we keep the data in order
	# dic['chromosome \t cpg position'] = coverage
	dic1 = defaultdict(lambda: 0)
	dic2 = defaultdict(lambda: 0)

	# Open and read
	with open(permeth1, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			try:
				name = line[0]
				meth = int(round(float(line[5]), 0))
				cov = int(line[4])
				if cov >= 50:
					dic1[name] = meth
			except:
				pass
		f.close()

	with open(permeth2, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			try:
				name = line[0]
				meth = int(round(float(line[5]), 0))
				cov = int(line[4])
				if cov >= 50:
					dic2[name] = meth
			except:
				pass
		f.close()


	# Put data into a numpy array 
	arrayoflists = []

		# Populate arrayoflists with [0,0,0,...]
		# Number of these arrays depend on unique cov of y (permeth2)
		# Number of 0's depends on unique cov of x (permeth1)

	for i in range(0, ybins+1):
		# This is done because we use log scale to view
		tmp = [1]*(ybins+1)
		arrayoflists.append(tmp)
	arrayoflists = np.array(arrayoflists)

	# Populate arrayoflists
	# 1st make set with all CpGs between both files
	totalcpgs = set(dic1.keys()) & set(dic2.keys())
	print len(totalcpgs)
	for pos in totalcpgs:
		cov1 = dic1[pos]
		cov2 = dic2[pos]
		if cov2 < ybins:
			if cov1 < ybins:
				y = ybins - cov2
				arrayoflists[y][cov1] += 1
			elif cov1 >= ybins:
				y = ybins - cov2
				arrayoflists[y][ybins] += 1
		
		elif cov2 >= ybins:
			if cov1 < ybins:
				arrayoflists[0][cov1] += 1

			elif cov1 >= ybins:
				arrayoflists[0][ybins] += 1

	# Also keep track of max values
	maxvar = 0
	for y in arrayoflists:
		for x in y:
			if x > maxvar:
				maxvar = x
	#print arrayoflists

	# Make labels:
	xlabels = []
	ylabels = []
	for i in range(0,ybins+1):
		xlabels.append(i)
	for i in range(0, ybins+1):
		ylabels.append(i)
	xlabels = map(str, xlabels)
	ylabels = sorted(ylabels, reverse=True)
	ylabels = map(str, ylabels)
	# Plot heatmap
	fig, ax = plt.subplots()
	#im = ax.imshow(arrayoflists)
	# Log norm
	im = ax.pcolor(arrayoflists[::-1], norm=LogNorm(vmin=1, vmax=maxvar))
	cbar = ax.figure.colorbar(im, ax=ax, **{})
	cbar.ax.set_ylabel('', rotation=-90, va="bottom")
	ax.set_xticks(map(lambda x: x+0.5, np.arange(len(xlabels))))
	ax.set_yticks(map(lambda y: y+0.5, np.arange(len(ylabels))))
	#ax.set_ylim(top=ybins-min(set(dic2.values()))+1)
	ax.set_xticklabels(xlabels)
	ax.set_yticklabels(ylabels[::-1], fontsize=2)
	# Rotate labels
	plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=2)
	# # in each box
	'''for i in range(len(ylabels)):
	    for j in range(len(xlabels)):
		text = ax.text(j, i, arrayoflists[i, j],
		               ha="center", va="center", color="w")
	'''
	ax.set_title("CpG Methylation")
	fig.tight_layout()
	plt.savefig(pdf, format='pdf')
	plt.close()
	pdf.close()
	
	

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='See comments in code')
        parser.add_argument('-x', dest='permeth1', type=str, help='Permeth file whose coverage to graph on the x axis')
        parser.add_argument('-y', dest='permeth2', type=str, help='Permeth file whose coverage to graph on the y axis')
        args=parser.parse_args()
        main(args)





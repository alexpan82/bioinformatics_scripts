# Takes 2 *longestTranscript.txt output files from klDiv_dat.py
# Plots the measure of choice b/t the files
from collections import defaultdict
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from sys import argv
from matplotlib import cm as cm
import matplotlib.colors as colors
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''Takes CLEAR *.dat files (make_dat.py) and returns genes that 
							differ significantly in mu-value across *.dat file(s)''')
	parser.add_argument('--x', required=True,
		        metavar='<x>',
		        help='File to plot on x-axis')
	parser.add_argument('--y', required=True,
		        metavar='<y>',
		        help='File to plot on y-axis')	   		       
	parser.add_argument('--output', required=True,
			metavar="<output>",
			help='File prefix for outputted png files')
	parser.add_argument('--measure', required=True,
			metavar="<measure>", default=0,
			help="Please specify 'kl', 'trunc', 'mu_diff', or 'exon_mu_diff'")
	
	args = parser.parse_args()

	x, y, outname, measure = args.x, args.y, args.output, args.measure
	
		
	# Combine tables and plot only shared genes
	df_list = []
	for f in [x, y]:
		df_list.append(pd.read_csv(f, sep="\t"))

	
	df_intersect = pd.concat(df_list, axis=1, join='inner')
	
	# Which measure to plot
	if measure == 'kl':
		xcol = 11
		axis_ranges = [0, 0.5, 0, 0.5]
		# Filter out genes that are do not have KL values in both files
		df_intersect = df_intersect.loc[(df_intersect.iloc[:,xcol] != -100) & (df_intersect.iloc[:,xcol+15] != -100)]
	elif measure == 'trunc':
		xcol = 12
		axis_ranges = [-1, 1, -1, 1]
	elif measure == 'mu_diff':
		xcol = 13
		axis_ranges = [-2, 2, -2, 2]
	elif measure == 'exon_mu_diff':
		xcol = 14
		axis_ranges = [-2, 2, -2, 2]
	else:
		print("Please specify 'kl', 'trunc', 'mu_diff', or 'exon_mu_diff'")
		quit()
	
	print('Intersect: %s' % df_intersect.shape[0])
	print('Truncated in %s: %s' % (x, df_intersect[df_intersect.iloc[:,7] < 0].shape[0]))
	print('Truncated in %s: %s' % (y, df_intersect[df_intersect.iloc[:,15] < 0].shape[0]))
	print('Truncated in both %s and %s: %s' % (x, y, df_intersect[(df_intersect.iloc[:,7] < 0) & (df_intersect.iloc[:,15] < 0)].shape[0]))
	df_union = pd.concat(df_list, axis=1, join='outer')
	print('Union: %s' % df_union.shape[0])

	# cols 7, 15 correspond to ratio_diff
	# Plot

	fig = plt.figure()
	plt.hexbin(df_intersect.iloc[:,xcol], df_intersect.iloc[:,xcol+15], cmap=cm.jet, bins='log', gridsize=1000)
	plt.axis(axis_ranges)
	cb = plt.colorbar()
	fig.savefig(outname + '.%s.png' % measure)


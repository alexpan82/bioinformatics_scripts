# Takes combined tsv from combExonDrop.py
# Also takes 2 column sufixes to compare
# Ex: --y HD170 --x MRD030
# Will compare exp_MRD030, exp_HD170, nonzero_ratio_MRD030, nonzero_ratio_HD170
# Combines tables and retains expression and exon drop out info for each sample

from collections import defaultdict
import pandas as pd
import argparse
import matplotlib.pyplot as plt

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''Takes CLEAR *.dat files (make_dat.py) and returns genes that 
							differ significantly in mu-value across *.dat file(s)''')
	parser.add_argument('--comb', required=True,
		        metavar='<comb>',
		        help='Combined tsv from combExonDrop.py')
	parser.add_argument('--y', required=True,
		        metavar='<y>',
		        help='A column SUFIX corresponding to one of the headers in the tsv file. Ex: nonzero_ratio_[SUFIX]')
	parser.add_argument('--x', required=True,
		        metavar='<x>',
		        help='A column SUFIX corresponding to one of the headers in the tsv file. Ex: nonzero_ratio_[SUFIX]')	   		       
	parser.add_argument('--output', required=True,
			metavar="<output>",
			help='Name of png file')
	parser.add_argument('--ymax', required=False, type=float, default=0,
			metavar="<ymax>",
			help='Limit on y-axis scale')
			
	args = parser.parse_args()
	print(args)
	
	fields = ['exons', 'length']
	for name in [args.y, args.x]:
		fields.append('tpm_' + name)
		fields.append('nonzero_ratio_' + name)

	df = pd.read_csv(args.comb, sep="\t", index_col=0)
	df = df[fields]
	
	# Filter out zeros first and calculate percentile ranks
	df = df[(df[fields[2]] > 0) | (df[fields[4]] > 0)]
	df['rank_' + args.y] = df[fields[2]].rank(pct=True)
	df['rank_' + args.x] = df[fields[4]].rank(pct=True)
	
	# Filter by exp/tpm
	df = df[(df[fields[2]] >= 1) & (df[fields[4]] >= 1)]
	# Filter by longest isoform (this is highly inefficient im sorry)
	isoform = defaultdict(lambda: 0)
	whichrow = defaultdict(lambda: None)
	for i in range(df.shape[0]):
		df_line = df.iloc[i]
		gene = df_line.name
		length = df_line['length']
		if isoform[gene] < length:
			isoform[gene] = length
			whichrow[gene] = i
	filter_rows = list(sorted(whichrow.values()))
	df = df.iloc[filter_rows]
	size = df.shape[0]

	# Plot hist
	df['ratio_diff'] = df[fields[3]] - df[fields[5]]
	hist = df['ratio_diff'].hist(bins=41, range=(-1,1))
	plt.title(args.output.split('.')[0] + ' (number of genes = %s)' % (str(size)))
	# set yaxis scale if desired
	if args.ymax != 0:
		plt.ylim([0, args.ymax])
	fig = hist.get_figure()
	fig.savefig(args.output)
	
	# Save as table
	df.to_csv(args.output.split('.')[0] + '.longestIso.txt', sep="\t")

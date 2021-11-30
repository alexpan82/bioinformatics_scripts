# KL divergence b/t each gene for 2 dat files
# Dat files from make_dat_printAll.py
# Ex: python klDiv_dat.py --p ex1.dat --q ex2.dat --outfile [prefix]

from collections import defaultdict
import pandas as pd
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def kl_div(p_dist, q_dist):
	# p_dist = np.array(p.split(','), dtype='float64')
	# q_dist = np.array(q.split(','), dtype='float64')
	
	p_dist = p_dist + 1
	q_dist = q_dist + 1
	sum_q_dist = sum(q_dist)
	
	p_dist = p_dist / sum(p_dist)
	q_dist = q_dist / sum(q_dist)
	
	kl = np.sum(np.where(p_dist != 0, p_dist * np.log(p_dist / q_dist), 0))
	return int(sum_q_dist), q_dist, kl
	
	
def simulate(sum_q, q_dist_norm, kl_limit, n):
	ratio = 0
	if n != 0:
		num_exon_bins = len(q_dist)
		store_kl = []
		for i in range(0, n):
			null_dist, edges = np.histogram(np.random.choice(num_exon_bins, size=sum_q, p=q_dist_norm), bins = num_exon_bins)
			null_dist = (null_dist + 1) / sum(null_dist)
			kl = np.sum(np.where(null_dist != 0, null_dist * np.log(null_dist / q_dist_norm), 0))
			store_kl.append(kl)
		ratio = len([i for i in store_kl if i >= kl_limit]) / n
	return ratio


def bin_cov(df, num_bins, utr):
	store_values = [[], [], [], [], [], [], [], []]
	for i in range(df.shape[0]):
		df_line = df.iloc[i]
		exon_cov = df_line['exon_cov']	# String with cov of each exon separated by ','
		num_exons = df_line['exons']
		num_nonzero = df_line['nonzero']
		mu_bp = df_line['mu']
		strand = bool(df_line['strand'])

		cov = np.array(exon_cov.split(','), dtype=int)	# Coverage dist as an array
		
		# If truncating utrs, need to remove 1st and last exons
		# Also need to recalculate num_exons and num_nonzero
		if utr is False:
			cov = cov[1:num_exons-1]
			num_exons -= 2
			num_nonzero = len([i for i in cov if i != 0])
		
		# cov = np.array_split(cov, num_bins)	# Split array into num_bins
		# cov = np.array(list(map(lambda x: sum(x), cov)), dtype='float64') # Elements are np.array_split are arrays that need to be consolidated
		norm_cov = cov
		
		# Normalize if possible
		sum_cov = sum(cov)
		if sum_cov > 0:
			norm_cov = np.array(list(map(lambda x: round(float(x) / sum_cov, 3), cov)), dtype='float64')
		# Truncate tailing zeros if there are less exons than specified bins
		#if num_exons < num_bins:
		#	norm_cov = norm_cov[0:num_exons]
		#	cov = cov[0:num_exons]
		
		consec_zeros = (None, None)
		fp_utr_percCov = None
		tp_utr_percCov = None
		# recalculate mu in terms of exon number instead of bp
		exonic_mu = 'None'
		if sum(cov) > 0:
			exonic_mu = ((2 / num_exons) * (sum(np.arange(num_exons) * cov) / sum(cov))) - 1
        # Calculate number of consecutive exons w/ no coverage on 5'- and 3'-ends
			consec_zeros = count_consecutive_zeros(cov)
			fp_utr_percCov = norm_cov[0]
			tp_utr_percCov = norm_cov[-1]
			if (strand is False):
				consec_zeros = (consec_zeros[1], consec_zeros[0],)
				fp_utr_percCov = norm_cov[-1]
				tp_utr_percCov = norm_cov[0]	
		store_values[0].append(cov)
		store_values[1].append(num_exons)
		store_values[2].append(num_nonzero)
		store_values[3].append(exonic_mu)
		store_values[4].append(consec_zeros[0])
		store_values[5].append(consec_zeros[1])
		store_values[6].append(fp_utr_percCov)
		store_values[7].append(tp_utr_percCov)

	# Replace values in dataframe
	df['exon_cov'] = store_values[0]
	df['exons'] = store_values[1]
	df['nonzero'] = store_values[2]
	df['exonic_mu'] = store_values[3]
	df['5p_consec_zero'] = store_values[4]
	df['3p_consec_zero'] = store_values[5]
	
	if utr is True:
		df['5p_utr_percCov'] = store_values[6]
		df['3p_utr_percCov'] = store_values[7]


def normalize(df):
	longest_iso = defaultdict(lambda: 0)
	exp_longest = defaultdict(lambda: 0)
	for i in range(df.shape[0]):
		df_line = df.iloc[i]
		gene = df_line['gene']
		length = df_line['length']
		rp = df_line['exp']     # reads/bp
		if longest_iso[gene] < length:
			longest_iso[gene] = length
			exp_longest[gene] = rp * 1000.  # reads per kilobase
	per_mil = sum(exp_longest.values()) / 1000000.  # per million scaling factor
	df['norm_exp'] = (df['exp'] * 1000.) / per_mil
	return None


def count_consecutive_zeros(cov_list):
	count_fprime, count_tprime = (0, 0)
	# Number of consecutive exons w/ zero coverage from 5'-end
	for c in cov_list:
		if c == 0:
			count_fprime += 1
		else:
			break
	# Number of consecutive exons w/ zero coverage from 3'-end
	for c in cov_list[::-1]:
		if c == 0:
			count_tprime += 1
		else:
			break
	return count_fprime, count_tprime
			
			
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''Takes CLEAR *.dat files (make_dat.py) and returns genes that 
							differ significantly in mu-value across *.dat file(s)''')
	parser.add_argument('--p', required=True,
		        metavar='<p>',
		        help='KL(P||Q): Dat file from make_dat_printAll.py')
	parser.add_argument('--q', required=True,
		        metavar='<q>',
		        help='KL(P||Q): Dat file from make_dat_printAll.py')	   		       
	parser.add_argument('--output', required=True,
			metavar="<output>",
			help='File prefix for of output txt and png files')
	parser.add_argument('--sim', required=False,
			metavar="<sim>", default=0, type=int,
			help='For each gene in Q, simulate randomly n times to determine significance of gene')
	parser.add_argument('--bin', required=False,
			metavar="<bin>", default=5, type=int,
			help='Limit on how many exons to consider for KL div test b/t 2 discrete distributions')
	parser.add_argument('--utr', required=False,
			metavar="<bin>", default=False, type=bool,
			help='By default, does NOT consider 1st and last exons. Will use *_noutr fields in dat file')
	parser.add_argument('--exp', required=False,
			metavar="<exp>", default=1, type=int,
			help='Performs KLD and plotting for only genes with expression >= the set value')
	
			
	args = parser.parse_args()
	print(args)
	
	# Set expression filter
	exp_filter = float(args.exp)
	
	if args.utr is True:
		fields = [0, 2, 4, 6, 7, 8, 10, 11]
	else:
		fields = [1, 3, 5, 6, 7, 9, 10, 11]
		
	names = ['exp', 'mu', 'length', 'gene', 'strand', 'nonzero', 'exons', 'exon_cov']
	df_list = []

	for dat_file in [args.p, args.q]:
		# Sort both files so each line is comparable
		# Fundamental assumption that dat files have the same number of lines (generated from same reference)
		df = pd.read_csv(dat_file, sep="\t", usecols=fields)
		df.columns = names
		df = df.sort_values(by=['gene', 'length', 'exons'])
		df['re_index'] = range(df.shape[0])
		df = df.rename(index=df['re_index'])
		
		# Regroup transcripts with more exons (bins) than allowed by args.bin
		# Removing UTRs is taken care of in bin_cov()
		bin_cov(df, args.bin, args.utr)
		normalize(df)
		
		# Add to list
		df_list.append(df)
		print(dat_file)
	
	# kl div for each gene/isoform that is expressed in both dat files
	kl_list = []
	kl_significance = []
	trunc_list = []
	longest_iso = defaultdict(lambda: 0)
	whichrow = defaultdict(lambda: None)
	for i in range(df_list[0].shape[0]):
		df_pline = df_list[0].iloc[i]
		df_qline = df_list[1].iloc[i]
		
		p_exons = df_pline['exon_cov']
		q_exons = df_qline['exon_cov']
		p_exp = df_pline['norm_exp']
		q_exp = df_qline['norm_exp']
		
		gene = df_pline['gene']
		length = df_pline['length']
		
		# Calc KL divergence if transcript >= exp_filter for both samples
		# if p_exp >= exp_filter and q_exp >= exp_filter:
		if p_exp >= exp_filter or q_exp >= exp_filter:
			sum_q_dist, q_dist, kl = kl_div(p_exons, q_exons)
			kl_list.append(kl)
			kl_significance.append(simulate(sum_q_dist, q_dist, kl, args.sim))
		else:
			kl_list.append(-100.)
			kl_significance.append(-100.)

		# Find rows with longest isoform
		if longest_iso[gene] < length:
			longest_iso[gene] = length
			whichrow[gene] = i
	
	# Write results to output
	adjust = 0
	out_fields = ['mu', 'exonic_mu', 'exp', 'norm_exp', 'nonzero', '5p_consec_zero', '3p_consec_zero']
	if args.utr is True:
		adjust = 2
		out_fields = ['mu', 'exonic_mu', 'exp', 'norm_exp', 'nonzero',
			'5p_consec_zero', '3p_consec_zero', '5p_utr_percCov', '3p_utr_percCov']
		
	result = pd.concat([df_list[0][['exons', 'length']],
				df_list[0][out_fields],
				df_list[1][out_fields]],
				axis=1, join='outer')

	result['kl'] = kl_list
	if args.sim != 0:
		result['simulate'] = kl_significance
	result = result.rename(index=df_list[0]['gene'])
	result = result.replace('None', np.NaN)
	#print(result)
	result['truncation'] = (result.iloc[:,6] - result.iloc[:,13+adjust]) / result.iloc[:,0]
	result['mu_diff'] = result.iloc[:,2].astype(float) - result.iloc[:,9+adjust].astype(float)
	result['exon_mu_diff'] = result.iloc[:,3] - result.iloc[:,10+adjust]
	result.to_csv(args.output + '.all.txt', sep="\t")
	
	# Filter by longest transcript for each gene
	filter_rows = list(sorted(whichrow.values()))
		
	# Plot filtered values
	kl_list = np.array(kl_list)
	kl_longest_iso = kl_list[filter_rows]	# Find value associated with longest isoform
	kl_longest_iso_clean = np.delete(kl_longest_iso, np.where(kl_longest_iso == -100))
	
	trunc_list = np.array(result['truncation'])
	trunc_longest_iso = trunc_list[filter_rows]
	trunc_longest_iso_clean = np.delete(trunc_longest_iso, np.where(kl_longest_iso == -100))
	
	p_exp_list = np.array(result.iloc[:,5])
	p_exp_list = p_exp_list[filter_rows]
	p_exp_list_clean = np.delete(p_exp_list, np.where(kl_longest_iso == -100))
	
	q_exp_list = np.array(result.iloc[:,12+adjust])
	q_exp_list = q_exp_list[filter_rows]
	q_exp_list_clean = np.delete(q_exp_list, np.where(kl_longest_iso == -100))
	
	result = result.iloc[filter_rows]
	result.to_csv(args.output + '.longestTranscript.txt', sep="\t")
	
	# Plot hist
	fig = plt.figure()
	plt.hist(kl_longest_iso_clean, bins=51, range=(0,10))
	plt.title('KL(%s || %s)' % (args.p, args.q))
	fig.savefig(args.output + '.kld.png')
	
	# Plot exon trun histogram
	'''
	fig = plt.figure()
	hist = result.loc[(result.iloc[:,2] >= exp_filter) & (result.iloc[:,5] >= exp_filter)]['truncation'].hist(bins=41, range=(-1,1))
	plt.title(args.output.split('.')[0])
	fig = hist.get_figure()
	fig.savefig(args.output + '.trunc.png')
	'''
	
	# Plot exon truncation vs KL divergence
	fig = plt.figure()
	plt.scatter(trunc_longest_iso_clean, kl_longest_iso_clean, alpha=0.3)
	plt.title('KL(%s || %s)' % (args.p, args.q))
	#plt.xlim([-1, 1])
	#plt.ylim([0, max(kl_longest_iso) + 1])
	fig.savefig(args.output + '.scatter.png')
	
	from matplotlib import cm as cm
	import matplotlib.colors as colors
	
	fig = plt.figure()
	plt.hexbin(trunc_longest_iso_clean, kl_longest_iso_clean, cmap=cm.jet, bins='log')
	plt.axis([-1, 1, 0, 8])
	plt.title('KL(%s || %s)' % (args.p, args.q))
	cb = plt.colorbar()
	#heatmap, xedges, yedges = np.histogram2d(trunc_longest_iso_clean, kl_longest_iso_clean, bins=50)
	# print(trunc_longest_iso_clean, kl_longest_iso_clean)
	# print(heatmap)
	#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	#plt.clf()
	#plt.imshow(heatmap.T, extent=extent, origin='lower')
	fig.savefig(args.output + '.heatmap.png')
	
	# Plot KL vs expression
	fig = plt.figure()
	
	num_pexp_bins = int(np.amax(p_exp_list_clean) / 4)
	num_kl_bins = int(np.amax(kl_longest_iso_clean)) * 10
	plt.hexbin(p_exp_list_clean, kl_longest_iso_clean, cmap=cm.jet, bins='log', gridsize=(num_pexp_bins, num_kl_bins))
	plt.axis([0, 300, 0, 10])
	plt.title('KL(%s || %s)' % (args.p, args.q))
	cb = plt.colorbar()
	fig.savefig(args.output + '.klvsP.png')
	
	fig = plt.figure()

	num_qexp_bins = int(np.amax(q_exp_list_clean) / 4)
	plt.hexbin(q_exp_list_clean, kl_longest_iso_clean, cmap=cm.jet, bins='log', gridsize=(num_qexp_bins, num_kl_bins))
	plt.axis([0, 300, 0, 10])
	plt.title('KL(%s || %s)' % (args.p, args.q))
	cb = plt.colorbar()
	fig.savefig(args.output + '.klvsQ.png')
	
	print('10th, 25th, 50th, 75th, 90th percentile of P: %s %s %s %s %s' % (np.percentile(p_exp_list_clean, 10),
	np.percentile(p_exp_list_clean, 25), np.percentile(p_exp_list_clean, 50), np.percentile(p_exp_list_clean, 75),
	np.percentile(p_exp_list_clean, 90)))
	
	print('10th, 25th, 50th, 75th, 90th percentile of Q: %s %s %s %s %s' % (np.percentile(q_exp_list_clean, 10),
	np.percentile(q_exp_list_clean, 25), np.percentile(q_exp_list_clean, 50), np.percentile(q_exp_list_clean, 75),
	np.percentile(q_exp_list_clean, 90)))
	
	'''
	heatmap, xedges, yedges = np.histogram2d(p_exp_list_clean, kl_longest_iso_clean, bins=1000)
	print(heatmap)
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	plt.clf()
	plt.imshow(heatmap.T, extent=extent, origin='lower', norm=colors.LogNorm(vmin=1, vmax=np.amax(heatmap)), aspect='equal')
	'''
	

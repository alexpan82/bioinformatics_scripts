# Takes CLEAR *.dat files (make_dat.py)
# Returns genes that differ significantly in mu-value across *.dat file(s)
# Chooses longest refseq isoform OR some other specified isoform

from collections import defaultdict
import pandas as pd
import argparse
import itertools
import matplotlib.pyplot as plt

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''Takes CLEAR *.dat files (make_dat.py) and returns genes that 
							differ significantly in mu-value across *.dat file(s)''')
	parser.add_argument('--dat', required=True,
		        metavar='<dat>', nargs='*',
		        help='CLEAR *.dat files (make_dat.py)')
		        
	parser.add_argument('--output', required=True,
			metavar="<output>",
			help='Name of output txt file')
			
	parser.add_argument('--filter', required=False,
			metavar="<filter>",
			help='By default chooses longest refseq isoform. If specified, will choose the refseq isoform specified in text file.')
			
			
	args = parser.parse_args()
	print(args)
	
	# dat_dict[file][gene] = [expression, mu, length, strand]
	dat_dict = defaultdict(lambda: defaultdict(lambda: [0, 'None', 0, 'None']))
	
	# 1) Read thru all dat files. Only choose longest isoform. Exclude genes with no expression
	# 2) Then turn into pandas df and calculate percentile of each filtered gene (based on exp)
	# 3) Split into 10 bins by expression. Calculate STD for each bin seperately
	df_dict = {}
	
	for dat_file in args.dat:
		with open(dat_file, 'r') as f:
			for line in f:
				line = line.strip().split('\t')
				exp = float(line[0])
				length = int(line[2])
				gene = line[3]
				strand = line[4]
				num_nonzero = int(line[5])
				num_exons = int(line[6])
				if (length > dat_dict[dat_file][gene][2]) and (exp > 0.0):
					mu = float(line[1])
					try:
						exon_bin = line[7]
					except:
						exon_bin = 'None'
					dat_dict[dat_file][gene] = [exp, mu, length, strand, num_nonzero, num_exons, exon_bin]

		# Filtered dict to pandas df
		df = pd.DataFrame.from_dict(dat_dict[dat_file], orient='index')
		df.columns = ['exp', 'mu', 'length', 'strand', 'num_nonzero', 'num_exons', 'exon_bin']
		df = df[df.mu != 'None']	# Do one last filtering just in case any null genes got throught
		df['pct_rank'] = df.exp.rank(pct = True)
		df['std'] = None	# Add empty STD column
		print(df.shape)
		
		# STD by bin
		for i in range(0, 50):
			lower_bound = i * 0.02
			upper_bound = (i+1) * 0.02
			tmp_df = df[(df.pct_rank > lower_bound) & (df.pct_rank <= upper_bound)]
			std = tmp_df.mu.std(axis=0)
			gene_names = tmp_df.index
			df.at[gene_names, 'std'] = std

		df_dict[dat_file] = df
	df_dict['L7365_WoyachJ_MRD030_BL_B_V1C_comb2.dat'].to_csv('hi.txt', sep="\t")
	
	# Intersect genes in each file
	set_list = []
	for df_dat in df_dict.keys():
		set_list.append(set(df_dict[df_dat].index))
	shared_genes = set.intersection(*set_list)
	del set_list
	print('Shared genes:', len(shared_genes))
	
	# Compare file pairs for shared genes
	# Only return genes whose mu values differ more than their collective std
	# Additionaly condition that the number of nonzero exons must change b/t samples
	#print('gene', 'mu_x', 'mu_y', 'std_x', 'std_y')\
	outdict = defaultdict(lambda: [])
	for gene in shared_genes:
		outdict['gene'].append(gene)
		for combo in itertools.combinations(args.dat[::-1], 2):
			mu_x = df_dict[combo[0]].loc[gene]['mu']
			mu_y = df_dict[combo[1]].loc[gene]['mu']
			std_x = df_dict[combo[0]].loc[gene]['std']
			std_y = df_dict[combo[1]].loc[gene]['std']
			exon_x = df_dict[combo[0]].loc[gene]['num_nonzero']
			exon_y = df_dict[combo[1]].loc[gene]['num_nonzero']
			exp_x = df_dict[combo[0]].loc[gene]['exp']
			exp_y = df_dict[combo[1]].loc[gene]['exp']
			num_exons = df_dict[combo[0]].loc[gene]['num_exons']
			
			mu_diff = abs(mu_x - mu_y)
			std_diff = abs(std_x + std_y)
			exon_diff = abs(exon_x - exon_y)
			
			if exp_x > 10 or exp_y > 10:
				if (exon_diff > 0):
					outdict[combo].append((round(mu_x,2), round(exp_x,0), round(mu_y,2), round(exp_y,0)))
					# print(gene, mu_x, mu_y, std_x, std_y)
				else:
					outdict[combo].append('NA')
			else:
				outdict[combo].append('NA')

			
	# Output to file		
	outdf = pd.DataFrame.from_dict(outdict)	
	outdf.to_csv(args.output, sep="\t")
	
	# Output mu std of each percentile of expression
	for i, f in enumerate(df_dict.keys()):
		if i == 0:
			ax = df_dict[f].sort_values(by=['pct_rank']).plot(x="pct_rank", y="std")
		else:
			df_dict[f].sort_values(by=['pct_rank']).plot(ax=ax, x="pct_rank", y="std")
		# fig = plot.get_figure()
	plt.gca().invert_xaxis()
	plt.xlabel('pct_rank')
	plt.ylabel("STD")
	#fig.savefig(f.split('.')[0] + ".std.png")
	fig = ax.get_figure()
	fig.savefig('test.png')

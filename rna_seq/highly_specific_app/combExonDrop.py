# Takes CLEAR *all.dat files (make_dat_printAll_v5.py)
# Combines tables and retains expression and exon drop out info for each sample

from collections import defaultdict
import pandas as pd
import argparse
import matplotlib.pyplot as plt

# Calculate tpm based on counts of longest isoforms
def tpm(df):
	longest_iso = defaultdict(lambda: 0)
	exp_longest = defaultdict(lambda: 0)
	for i in range(df.shape[0]):
		df_line = df.iloc[i]
		gene = df_line['gene']
		length = df_line['length']
		rp = df_line['exp']	# reads/bp
		if longest_iso[gene] < length:
			longest_iso[gene] = length
			exp_longest[gene] = rp * 1000.	# reads per kilobase
	per_mil = sum(exp_longest.values()) / 1000000.	# per million scaling factor
	df['tpm'] = (df['exp'] * 1000.) / per_mil
	return None

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''Takes CLEAR *.dat files (make_dat.py) and returns genes that 
							differ significantly in mu-value across *.dat file(s)''')
	parser.add_argument('--dat', required=True,
		        metavar='<dat>', nargs='*',
		        help='CLEAR *.dat files (make_dat_printAll_v5.py)')
		        
	parser.add_argument('--output', required=True,
			metavar="<output>",
			help='Name of output txt file')
	
	parser.add_argument('--utr', required=True,
			metavar="<output>",
			help='T to include UTRs, F to remove UTRs. Will use *_noutr fields in dat file')
			
			
	args = parser.parse_args()
	print(args)
	
	# fields = [exp, length, gene, non_zero_exons, exons]
	if args.utr == 'T':
		fields = [0, 4, 6, 8, 10]
	elif args.utr == 'F':
		fields = [1, 5, 6, 9, 10]
	else:
		print('--utr T to include UTRs, F to remove UTRs')
		quit()

	names = ['exp', 'length', 'gene', 'non_zero_exons', 'exons']
	df_list = [None, None, None]

	for dat_file in args.dat:
		df = pd.read_csv(dat_file, sep="\t", usecols=fields)
		df.columns = names
		
		# If removing UTRs, remove transcripts with < 2 exons
		if args.utr == 'F':
			df['exons'] = df['exons'] - 2
			df = df[(df['exons'] > 0)]
			
		df['nonzero_ratio'] = df['non_zero_exons'] / df['exons']
		df = df.sort_values(by=['gene', 'length'])
		df['re_index'] = range(df.shape[0])
		# index_names = df['gene'] + '_' + df['length'].astype('str')
		df = df.rename(index=df['re_index'])
		# tpm
		tpm(df)
		
		df_list[0] = df['gene']
		df_list[1] = df['exons']
		df_list[2] = df['length']
		
		
		df = df[['tpm', 'nonzero_ratio']]
		dat_name = dat_file.split('.')[0]
		df = df.rename(columns={'tpm': 'tpm_' + dat_name, 'nonzero_ratio': 'nonzero_ratio_' + dat_name})
		df_list.append(df)
	
	result = pd.concat(df_list[1:], axis=1, join='outer')
	result = result.rename(index=df_list[0])
	result.to_csv(args.output, sep="\t")

